if(!require("choisycutoff")){
  devtools::install_github("shengqh/cutoff")
}
library("choisycutoff")
library("MASS")
library("reshape2")
library("ggplot2")
library("Seurat")
library("gridExtra")
library("ggExtra")

# ----

is_windows<-function(){
  return(.Platform$OS.type == "windows")
}

my_startval <- function(values,D1="normal",D2="normal",cutoff_point=0) {
  if(cutoff_point >0){
    thresh = cutoff_point
  }else{
    den <- tryCatch(
      expr = {
        density(values, bw="SJ")
      },
      error = function(e){ 
        density(values)
      }
    )
    w=1
    x=den$x
    y=den$y
    y.smooth=den$y
    n <- length(y.smooth)
    y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
    delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
    i.max <- which(delta <= 0) + w
    res=data.frame(x=x[i.max], i=i.max, y=y[i.max])
    res=res[res$x>0,]
    res=res[order(res$y, decreasing = T),]
    
    if(nrow(res)>2){
      res=res[1:2,]
    }
    
    xx=x[x>min(res$x) & x<max(res$x)]
    yy=y[x>min(res$x) & x<max(res$x)]
    yy.min=min(yy)
    ii.min=which(yy==yy.min)
    thresh=xx[ii.min]
  }
  
  sel <- values<thresh
  data1 <- values[sel]
  data2 <- values[!sel]
  lambda <- length(data1)/length(values)
  param1 <- MASS::fitdistr(data1,D1)$est
  param2 <- MASS::fitdistr(data2,D2)$est
  out <- c(param1,param2,lambda)
  names(out) <- c("mu1","sigma1","mu2","sigma2","lambda")
  return(out)
}

# ----

my_em<-function(values, data_name="em", D1="normal", D2="normal", t=1e-64, cutoff_point=0, max_iteration=1000){
  start <- as.list(my_startval(values, D1, D2, cutoff_point))
  
  D1b <- choisycutoff:::hash[[D1]]
  D2b <- choisycutoff:::hash[[D2]]
  lambda0 <- 0
  with(start, {
    iteration = 0
    while (abs(lambda0 - mean(lambda)) > t) {
      lambda <- mean(lambda)
      lambda0 <- lambda
      distr1 <- lambda * D1b(values, mu1, sigma1)
      distr2 <- (1 - lambda) * D2b(values, mu2, sigma2)
      lambda <- distr1/(distr1 + distr2)
      mLL2 <- function(mu1, sigma1, mu2, sigma2) return(choisycutoff:::mLL(mu1, 
                                                                           sigma1, mu2, sigma2, lambda, values, D1b, D2b))
      start <- as.list(log(c(mu1 = mu1, sigma1 = sigma1, 
                             mu2 = mu2, sigma2 = sigma2)))
      out <- bbmle::mle2(mLL2, start, "Nelder-Mead")
      coef <- out@coef
      coef_n <- names(coef)
      names(coef) <- NULL
      for (i in 1:4) assign(coef_n[i], exp(coef[i]))
      iteration = iteration + 1
      if(iteration >= max_iteration){
        warning(paste0("reach ", max_iteration, " iterations, break."))
        break
      }
    }
    out <- list(lambda = lambda, param = exp(out@coef), D1 = D1, 
                D2 = D2, deviance = out@min, data = values, data_name = data_name, 
                out = out, t = t)
    class(out) <- "em"
    return(out)
  })
}

# ----

my_cutoff<-function (my_out, t = 1e-64, nb = 10, distr = 2, type1 = 0.05, level = 0.95) 
{
  coef <- my_out$out@coef
  the_names <- names(coef)
  coef <- exp(mc2d::rmultinormal(nb, coef, as.vector(my_out$out@vcov)))
  coef <- as.list(data.frame(t(coef)))
  coef <- lapply(coef, function(x) {
    names(x) <- the_names
    return(as.list(x))
  })
  out <- sapply(coef, function(x) choisycutoff:::lci0(x, mean(my_out$lambda), 
                                                      choisycutoff:::hash[[my_out$D1]], choisycutoff:::hash[[my_out$D2]], my_out$data, t))
  lambda <- rnorm(nb, out[1, ], out[2, ])
  coef <- sapply(coef, function(x) unlist(x))
  the_names <- c(rownames(coef), "lambda")
  coef <- rbind(coef, lambda)
  coef <- lapply(as.data.frame(coef), function(x) {
    names(x) <- the_names
    return(as.list(x))
  })
  out <- sapply(coef, function(x) with(x, choisycutoff:::cutoff0(mu1, 
                                                                 sigma1, mu2, sigma2, lambda, my_out$D1, my_out$D2, distr, type1)))
  out <- MASS::fitdistr(out, "normal")
  the_mean <- out$estimate["mean"]
  level <- (1 - level)/2
  level <- c(level, 1 - level)
  ci <- the_mean + qt(level, Inf) * out$sd["mean"]
  out <- c(the_mean, ci)
  names(out) <- c("Estimate", paste(100 * level, "%"))
  return(out)
}

# ----

get_cutoff<-function(values, prefix=NULL, cutoff_startval=0){
  my_out <- my_em(values,"normal","normal", cutoff_point=cutoff_startval)
  
  if(!is.null(prefix)){
    saveRDS(my_out, paste0(prefix, ".em.rds"))
  }

  cut_off <- tryCatch({
      my_cutoff(my_out)
    }, error=function(e){
      if(cutoff_startval == 0){
        stop(paste0("failed to find cutoff due to: ", e))
      }else{
        print(paste0("failed to find cutoff due to: ", e, ", use start value ", cutoff_startval, " as cutoff"))
        return(cutoff_startval)
      }
    }
  )

  if(!is.null(prefix)){
    png(paste0(prefix, ".cutoff.png"), width=2000, height=1600, res=300)
    hist(values,200,F,xlab="concentration",ylab="density", main=NULL,col="grey")
    lines(density(values),lwd=1.5,col="blue")
    lines(my_out,lwd=1.5,col="red")
    if(cutoff_startval != 0){
      abline(v=cutoff_startval,lwd=1.5,col="green")
    }
    abline(v=cut_off[1],lwd=1.5,col="brown")
    dev.off()
  }
  
  return(cut_off[1])
}

rplot<-function(object, features, assay, identName, withAllCells=FALSE, n_row=1){
  DefaultAssay(object = object) <- assay
  data <- FetchData(object = object, vars = c(features, identName))
  mdata<-melt(data, id.vars=identName)
  if (withAllCells) {
    mdata2<-mdata
    mdata2[,1] = "All cells"
    mdata<-rbind(mdata, mdata2)
  }
  
  gfinal=list()
  for(feature in features){
    ddata=mdata[mdata$variable==feature,]
    g<-ggplot(ddata, aes_string(x="value")) + 
      geom_histogram(aes(y=..density..), bins=50, colour="black", fill="white", position="identity") + 
      geom_density(color="red") +
      xlab(feature) + theme_bw()+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=24),
            axis.title.y=element_blank())
    if(length(unique(mdata[,identName])) > 1){
      g<-g+facet_grid(reformulate(".", identName), scale="free_y") + 
        theme(strip.background=element_rect(colour="black", fill=NA),
              strip.text = element_text(size = 24))
    }
    gfinal = append(gfinal, list(g))
  }
  grid.arrange(grobs=gfinal, nrow=n_row)
}


# ----

get_object<-function(counts) {
  obj <- CreateSeuratObject(counts = counts, assay="HTO")
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
  DefaultAssay(object = obj) <- "HTO"
  
  return(obj)
}

read_hto<-function(rdsfile) {
  htos<-readRDS(rdsfile)
  return(get_object(htos))
}

output_post_classification<-function(obj, output_prefix, demultiplex_column="scDemultiplex"){
  tagnames=rownames(obj[["HTO"]])
  
  hto_names=unique(obj@meta.data[,demultiplex_column])
  a_hto_names=hto_names[!(hto_names %in% c("Doublet","Negative"))]
  a_hto_names=a_hto_names[order(a_hto_names)]
  hto_names=c(a_hto_names, "Negative", "Doublet")
  obj$HTO_classification=factor(obj$HTO_classification, levels=hto_names)
  
  width=max(1600, length(tagnames) * 1000)
  
  Idents(obj) <- "HTO_classification"
  png(paste0(output_prefix, ".class.ridge.png"), width=width, height=max(1400, length(tagnames) * 300), res=300)
  print(RidgePlot(obj, assay = "HTO", features = tagnames, ncol = length(tagnames)))
  dev.off()
  
  png(paste0(output_prefix, ".class.dist.png"), width=width, height=max(1400, length(tagnames) * 500), res=300)
  rplot(obj, assay = "HTO", features = tagnames, identName="HTO_classification")
  dev.off()
  
  if (length(tagnames) == 2) {
    png(paste0(output_prefix, ".class.point.png"), width=2000, height=1800, res=300)
    print(FeatureScatter(object = obj, feature1 = tagnames[1], feature2 = tagnames[2],group.by="HTO_classification"))
    dev.off()
  }
  
  tmat=data.frame(t(data.frame(obj@assays$HTO@counts)))
  rownames(tmat)=colnames(obj)
  tmat$HTO = unlist(obj$HTO_classification)
  tmat$HTO.global = unlist(obj$HTO_classification.global)
  write.csv(tmat, file=paste0(output_prefix, ".csv"))
  
  if(length(tagnames) >= 2) {
    obj<-ScaleData(obj, features=tagnames, verbose=FALSE)
    obj<-RunPCA(obj, features=tagnames, approx=FALSE)    
    obj<-RunUMAP(obj, features=tagnames, slot="scale.data")
    
    png(paste0(output_prefix, ".umap.class.png"), width=2000, height=1800, res=300)
    g<-DimPlot(obj, reduction = "umap", group.by="HTO_classification")
    print(g)
    dev.off()
    
    nwidth=ceiling(sqrt(length(tagnames)))
    nheight=ceiling(length(tagnames)/nwidth)
    png(paste0(output_prefix, ".umap.tag.png"), width=nwidth*1500, height=1500*nheight, res=300)
    g<-FeaturePlot(obj, features=tagnames, reduction = "umap")
    print(g)
    dev.off()
    
    cols=rep("gray", length(hto_names))
    names(cols)=hto_names
    cols[['Negative']]="blue"
    cols[["Doublet"]]="red"
    
    png(paste0(output_prefix, ".umap.all.png"), width=2000, height=1800, res=300)
    g<-DimPlot(obj, reduction = "umap", label=T, group.by="HTO_classification", order=c("Negative", "Doublet"))+
      scale_color_manual(values=cols)
    print(g)
    dev.off()
  }
}


# ---


build_summary<-function(allres, output_prefix, nameMapFile=NA){
  colnames(allres)<-c("File", "Sample")
  dat=lapply(allres$File, function(x){
    dd=read.csv(x, check.names=F)
    table(dd$HTO.global)
  })
  dat.all=do.call(rbind, dat)
  rownames(dat.all)=allres$Sample
  write.csv(dat.all, file=paste0(output_prefix, ".summary.csv"), quote=F)
  
  mdat=reshape2::melt(dat.all)
  colnames(mdat)=c("Sample", "Class", "Cell")
  
  png(paste0(output_prefix, ".summary.png"), width=1600, height=1200, res=300)
  g<-ggplot(mdat, aes(x=Sample, y=Cell, fill=Class, label=Cell)) + geom_bar(position="stack", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(g)
  dev.off()
}




