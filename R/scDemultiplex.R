library(Seurat)
library(ggplot2)
library(zoo)
library(reshape2)
library(gridExtra)
library(ggExtra)
library(TailRank) # dbb # Beta-Binomial Distribution
library(edgeR)
library(dirmult)
library(MGLM) # ddirmn
library(parallel)
library(tictoc)

check_mc_cores<-function(mc.cores) {  
  if(.Platform$OS.type == "windows") {
    mc.cores=1
  }else{
    mc.cores=min(parallel::detectCores() - 1, max(1, mc.cores))
  }
  return(mc.cores)
}

do_cutoff<-function(tagname, data, output_prefix, cutoff_startval){
  values=data[,tagname]
  values=values[values>0] # remove count = 0

  if(is.list(cutoff_startval)){
    if(tagname %in% names(cutoff_startval)){
      cur_cutoff = cutoff_startval[tagname]
    }else{
      cur_cutoff = 0
    }
  }else{
    cur_cutoff = cutoff_startval
  }
  cat(paste0("get cutoff of ", tagname, " ...\n"))
  cutoff=get_cutoff(values, paste0(output_prefix, "_", tagname), cur_cutoff)
  return(cutoff)
}

#' @export
demulti_cutoff<-function(counts, output_prefix, cutoff_startval=0, mc.cores=1, cutoff_list=NULL){
  if(is(counts,"Seurat")){
    obj=counts
  }else{
    obj<-get_object(counts)
  }
  
  output_file = paste0(output_prefix, ".csv")
  tagnames = rownames(obj[["HTO"]])
  tagnames = tagnames[order(tagnames)]
  
  data <- FetchData(object=obj, vars=tagnames)
  
  if(!is.null(cutoff_list)){
    if(!is.list(cutoff_list)){
      stop(paste0("cutoff_list has to be named list which conatins all tagnames: ", paste0(tagnames, collapse = ",")))
    }
    if(!all(tagnames %in% names(cutoff_list))){
      stop(paste0("cutoff_list has to be named list which conatins all tagnames: ", paste0(tagnames, collapse = ",")))
    }
  }else{
    if(!require("choisycutoff")){
      BiocManager::install('shengqh/cutoff')
      library(choisycutoff)
    }
    mc.cores<-check_mc_cores(mc.cores)
    cutoff_list<-unlist(parallel::mclapply(tagnames, do_cutoff, data = data, output_prefix = output_prefix, cutoff_startval = cutoff_startval, mc.cores=mc.cores))
    names(cutoff_list) = tagnames
  }
  saveRDS(cutoff_list, paste0(output_prefix, ".cutoff_list.rds"))

  print(cutoff_list)

  tagname=tagnames[1]  
  for (tagname in tagnames) {
    cutoff = cutoff_list[tagname]
    data[,paste0(tagname,"_pos")] = ifelse(data[,tagname]>cutoff, tagname, "Negative")
  }
  
  cat(paste0("get classification ...\n"))
  class_names=paste0(tagnames, "_pos")
  data$scDemultiplex_cutoff=unlist(apply(data, 1, function(x){
    xx=unique(x[class_names])
    if (length(xx) > 1){
      xx = xx[xx != "Negative"]
      if (length(xx) > 1) {
        return("Doublet")
      }
    }
    return(xx)
  }))
  
  obj$scDemultiplex_cutoff = factor(data$scDemultiplex_cutoff, levels=c(tagnames, "Negative", "Doublet"))
  obj$scDemultiplex_cutoff.global = ifelse(obj$scDemultiplex_cutoff %in% c("Negative", "Doublet"), as.character(obj$scDemultiplex_cutoff), "Singlet")
  return(obj)
}
  
#Dirichlet-Multinomial Distribution
estimate_alpha<-function(name, taglist){
  x <-taglist[[name]]
  p.tu <- goodTuringProportions(colSums(x))
  print(paste0("    dirmult ", name, " ..."))
  theta <- dirmult(x, trace=F)$theta
  alpha_t <- (1-theta)/theta
  alpha_est <- alpha_t*p.tu
  return(alpha_est)
}

#' @export
demulti_refine<-function(obj, output_prefix="demulti_refine", p.cut=0.001, iterations=10, init_column="scDemultiplex_cutoff", mc.cores=1, refine_negative_doublet_only=FALSE){
  mc.cores<-check_mc_cores(mc.cores)

  dd <- obj[["HTO"]]@counts
  dd <- t(as.matrix(dd)) # 8193   12
  dd <- as.data.frame(dd)
  
  tag.var <- names(dd)
  
  nn.tag <- length(tag.var)
  dd$total.count <- apply(dd, 1, sum)
  NN <- dd$total.count
  
  start_time2 <- Sys.time()
  
  dd$HTO_classification <- unlist(obj[[init_column]])

  df<-data.frame(table(dd$HTO_classification))
  df$Iteration=0
  
  #cl <- makeCluster(nn.tag)  
  #registerDoParallel(cl) 
  
  kk=1
  for(kk in 1:iterations){
    print(paste0("refine iteration ", kk))

    lastClassification = dd$HTO_classification
    
    taglist <- split(dd[tag.var], dd$HTO_classification)
    taglist <- taglist[! names(taglist) %in% c("Negative","Doublet")]
    
    print("  estimate alpha ...")
    tic()
    out.alpha.est <- parallel::mclapply(names(taglist), estimate_alpha, taglist=taglist, mc.cores=mc.cores)
    toc()
    names(out.alpha.est)<-names(taglist)
    
    # ----
    print("  calculate pvalue ...")
    for(j in 1:nn.tag){
      alpha.est <- out.alpha.est[[tag.var[j]]]
      uu <- alpha.est[tag.var[j],]
      vv <- sum(alpha.est) - alpha.est[tag.var[j],]
      
      var.new <- paste(tag.var[j], ".pbb", sep = "")
      dd[,var.new] <- NA
      dd[,var.new] <- sapply(1:nrow(dd), function(x) pbb(dd[x, tag.var[j]], NN[x], uu, vv))
      dd[,var.new] <- ifelse(dd[,var.new] < 0.5, dd[,var.new], 1 - dd[,var.new])
      
      var.new1 <- paste(var.new, ".ps.pvalue.fdr", sep = "")
      dd[,var.new1] <- p.adjust(dd[,var.new], method = "BH")
      
      rm(var.new); rm(var.new1)
    }
    rm(j)
    
    # ----
    
    print("  assign category ...")
    dd$HTO_classification.comb2 <- NA
    dd$HTO_classification.list2 <- NA
    #p.cut <- 0.025
    #varp <- paste(tag.var, ".pbb", sep = "")
    #p.cut <- 0.001
    varp <- paste(tag.var, ".pbb.ps.pvalue.fdr", sep = "")
    
    for(i in 1:nrow(dd)){
      out <- dd[i,varp]
      out2 <- tag.var[which(out > p.cut)]
      dd$HTO_classification.list2[i] <- length(out2)
      dd$HTO_classification.comb2[i] <- paste(out2, collapse = "; ")
      rm(out); rm(out2)
    }
    rm(i)
    
    dd$HTO_classification.comb2[which(dd$HTO_classification.list2 == 0)] <- "Negative"
    
    # ----
    
    for(i in 1:length(tag.var)){
      if(refine_negative_doublet_only){
        dd$HTO_classification[which(dd$HTO_classification.list2 == 1 & 
                                    dd$HTO_classification %in% c("Negative", "Doublet") &
                                    dd$HTO_classification.comb2 == tag.var[i])] <- tag.var[i]
      }else{
        dd$HTO_classification[which(dd$HTO_classification.list2 == 1 & 
                                    dd$HTO_classification.comb2 == tag.var[i])] <- tag.var[i]
      }
    }
    rm(i)

    if(all(lastClassification == dd$HTO_classification)){
      break
    }

    cur_df<-data.frame(table(dd$HTO_classification))
    cur_df$Iteration=kk   

    df<-rbind(df, cur_df) 
  }
  rm(kk)

  mdf<-dcast(df, Var1~Iteration, value.var="Freq")
  write.csv(mdf, paste0(output_prefix, ".iteration.csv"))
  
  end_time2 <- Sys.time()
  time_run2 <- end_time2 - start_time2
  time_run2 # Time difference of 23.76274 mins
  
  # ----
  
  obj$scDemultiplex = factor(dd$HTO_classification, levels=c(tag.var, "Negative", "Doublet"))
  obj$scDemultiplex.global = ifelse(obj$scDemultiplex %in% c("Negative", "Doublet"), as.character(obj$scDemultiplex), "Singlet")
  
  return(obj)
}

#' @export
hto_umap<-function(obj){
  tagnames<-rownames(obj)
  DefaultAssay(obj)<-"HTO"
  obj <- ScaleData(obj, features=tagnames)
  obj <- RunPCA(obj, features=tagnames, approx = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:length(tagnames))
  return(obj)
}

#' @export
hto_plot<-function(obj, output_prefix, group.by){
  tagnames<-rownames(obj)

  g<-DimPlot(obj, reduction = "umap", group.by = group.by)
  png(paste0(output_prefix, ".umap.class.png"), width=2000, height=1800, res=300)
  print(g)
  dev.off()
  
  width=max(1600, length(tagnames) * 1000)
  
  Idents(obj) <- group.by
  png(paste0(output_prefix, ".class.ridge.png"), width=width, height=max(1400, length(tagnames) * 300), res=300)
  print(RidgePlot(obj, assay = "HTO", features = tagnames, ncol = length(tagnames)))
  dev.off()
  
  png(paste0(output_prefix, ".class.dist.png"), width=width, height=max(1400, length(tagnames) * 500), res=300)
  rplot(obj, assay = "HTO", features = tagnames, identName=group.by)
  dev.off()
  
  if (length(tagnames) == 2) {
    png(paste0(output_prefix, ".class.point.png"), width=2000, height=1800, res=300)
    print(FeatureScatter(object = obj, feature1 = tagnames[1], feature2 = tagnames[2],group.by=group.by))
    dev.off()
  }
  
  nwidth=ceiling(sqrt(length(tagnames)))
  nheight=ceiling(length(tagnames)/nwidth)
  png(paste0(output_prefix, ".umap.tag.png"), width=nwidth*1500, height=1500*nheight, res=300)
  g<-FeaturePlot(obj, features=tagnames, reduction = "umap")
  print(g)
  dev.off()
  
  cols=rep("gray", length(tagnames))
  names(cols)=tagnames
  cols[['Negative']]="blue"
  cols[["Doublet"]]="red"
  
  png(paste0(output_prefix, ".umap.nd.png"), width=2000, height=1800, res=300)
  g<-DimPlot(obj, reduction = "umap", label=T, group.by=group.by, order=c("Negative", "Doublet"))+
    scale_color_manual(values=cols)
  print(g)
  dev.off()
  
# 
#   hto_names <- c(tagnames, "Singlet", "Doublet", "Doublet to Singlet", "Negative to Singlet")
#   cols=rep("gray", length(hto_names))
#   names(cols)=hto_names
#   cols[['Negative']]="blue"
#   cols[["Doublet"]]="red"
#   cols[["Doublet to Singlet"]]="green"
#   cols[["Negative to Singlet"]]="gold"
#   
#   png(paste0(output_prefix, ".umap.all.png"), width=2000, height=1800, res=300)
#   g<-DimPlot(obj, reduction = "umap", label=T, group.by="HTO_classification3", order=c("Negative", "Doublet", "Doublet to Singlet", "Negative to Singlet"))+
#     scale_color_manual(values=cols)
#   print(g)
#   dev.off()
  
  return(obj)
}

#' @export
hto_findFoldChange<-function(obj, col) {
  Idents(obj)<-col
  markers<-FindAllMarkers(obj, pseudocount.use = FALSE,only.pos = TRUE)
  markers<-markers[!(markers$cluster %in% c("Negative", "Doublet")),]
  markers<-markers[order(markers$gene),c("gene", "avg_log2FC"), drop=F]
  rownames(markers)<-markers$gene
  markers<-markers[,"avg_log2FC", drop=F]
  colnames(markers)=col
  return(markers)
}

#' @export
hto_findMultiFoldChanges<-function(obj, cols) {
  alltb<-NULL

  for (col in cols){
    markers<-scDemultiplex_findFoldChange(obj, col)
    if(is.null(alltb)){
      alltb<-markers
    }else{
      alltb<-cbind(alltb, markers)
    }
  }
  
  return(alltb)
}