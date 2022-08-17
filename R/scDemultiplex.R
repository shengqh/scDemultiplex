library(Seurat)
library(ggplot2)
library(choisycutoff)
library(zoo)
library(reshape2)
library(gridExtra)
library(ggExtra)
library(TailRank) # dbb # Beta-Binomial Distribution
library(edgeR)
library(dirmult)
library(MGLM) # ddirmn

scDemultiplex_by_cutoff<-function(counts, output_prefix, cutoff_startval=0){
  obj<-get_object(counts)
  
  output_file = paste0(output_prefix, ".csv")
  tagnames = rownames(obj[["HTO"]])
  
  data <- FetchData(object=obj, vars=tagnames)
  
  tagname=tagnames[1]  
  for (tagname in tagnames) {
    values=data[,tagname]
    values=values[values>0] # remove count = 0
    cat(paste0("get cutoff of ", tagname, " ...\n"))
    cutoff=get_cutoff(values, paste0(output_prefix, "_", tagname), cutoff_startval)
    data[,paste0(tagname,"_pos")] = ifelse(data[,tagname]>cutoff, tagname, "Negative")
  }
  
  
  cat(paste0("get classification ...\n"))
  class_names=paste0(tagnames, "_pos")
  data$HTO_classification=unlist(apply(data, 1, function(x){
    xx=unique(x[class_names])
    if (length(xx) > 1){
      xx = xx[xx != "Negative"]
      if (length(xx) > 1) {
        return("Doublet")
      }
    }
    return(xx)
  }))
  
  
  data$HTO_classification.global=unlist(apply(data, 1, function(x){
    xx=unique(x[class_names])
    if (length(xx) > 1){
      xx = xx[xx != "Negative"]
      if (length(xx) > 1) {
        return("Doublet")
      }else{
        return("Singlet")
      }
    }
    return("Negative")
  }))
  
  
  obj$HTO_classification = data$HTO_classification
  obj$HTO_classification2 = factor(obj$HTO_classification, levels=c(tagnames, "Negative", "Doublet"))
  obj$HTO_classification.global = data$HTO_classification.global
  return(obj)
}

scDemultiplex_by_refine<-function(obj, p.cut=0.0001){
  dd <- obj[["HTO"]]@counts
  dd <- t(as.matrix(dd)) # 8533   12
  dd <- as.data.frame(dd)
  
  tag.var <- names(dd)
  tag.var1 <- c(tag.var, "Negative")
  
  dd$total.count <- apply(dd, 1, sum)
  
  dd$id <- row.names(dd)
  
  # ----
  
  start_time2 <- Sys.time()
  
  dd$HTO_classification.tiger <- dd$HTO_classification.r0 <- dd$HTO_classification <- obj$HTO_classification
  dd$HTO_classification.global.tiger <- dd$HTO_classification.global.r0 <- dd$HTO_classification.global <- obj$HTO_classification.global
  
  for(kk in 1:10){
    print(paste0("iteration ", kk))
    # start_time1 <- Sys.time()
    
    out.alpha.est <- matrix(NA, 2, length(tag.var))
    for(i in 1:length(tag.var)){
      
      dds <- subset(dd, HTO_classification == tag.var[i], select = tag.var) # HTO_classification
      tag.sum <- apply(dds[,tag.var], 2, sum)
      tag.sum <- c(tag.sum[i], sum(tag.sum[-i])) # make it to k=2
      p.tu <- goodTuringProportions(tag.sum)
      dds2 <- data.frame(x1=dds[,tag.var[i]],
                         x2=apply(dds[,tag.var[-i]], 1, sum))
      fit <- dirmult(dds2, init = p.tu, epsilon=10^(-4), trace = FALSE)
      theta.est <- fit$theta
      alpha.est <- p.tu*(1-theta.est)/theta.est
      out.alpha.est[, i] <- alpha.est
      
      rm(dds); rm(tag.sum); rm(p.tu); rm(fit); rm(theta.est)
    }
    rm(i)
    
    # end_time1 <- Sys.time()
    # time_run1 <- end_time1 - start_time1
    # time_run1 # 
    
    # ----
    
    # start_time2 <- Sys.time()
    
    NN <- dd$total.count
    nn.tag <- length(tag.var) 
    for(j in 1:nn.tag){
      
      alpha.est <- out.alpha.est[,j]
      uu <- alpha.est[1]
      vv <- alpha.est[2]
      
      var.new <- paste(tag.var[j], ".pbb", sep = "")
      dd[,var.new] <- NA
      dd[,var.new] <- sapply(1:nrow(dd), function(x) pbb(dd[x, tag.var[j]], NN[x], uu, vv))
      dd[,var.new] <- ifelse(dd[,var.new] < 0.5, dd[,var.new], 1 - dd[,var.new])
      
      var.run <- paste(tag.var[j], ".pbb.r", kk, sep = "")
      dd[,var.run] <- dd[,var.new] 
      rm(var.run)
      
      var.new1 <- paste(var.new, ".ps.pvalue.fdr", sep = "")
      dd[,var.new1] <- p.adjust(dd[,var.new], method = "BH")
      
      var.run <- paste(var.new1, ".r", kk, sep = "")
      dd[,var.run] <- dd[,var.new1]
      rm(var.run)
      
      rm(var.new); rm(var.new1)
    }
    rm(j)
    
    # end_time2 <- Sys.time()
    # time_run2 <- end_time2 - start_time2
    # time_run2 # 
    
    # ----
    
    dd$HTO_classification.comb2 <- NA
    dd$HTO_classification.list2 <- NA
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
    
    var.hot <- paste("HTO_classification.comb2.r", kk, sep ="")
    dd[,var.hot] <- dd$HTO_classification.comb2
    rm(var.hot)
    
    var.hot <- paste("HTO_classification.list2.r", kk, sep ="")
    dd[,var.hot] <- dd$HTO_classification.list2
    rm(var.hot)
    
    # ----
    
    for(i in 1:length(tag.var)){
      dd$HTO_classification[which(dd$HTO_classification.list2 == 1 & 
                                    dd$HTO_classification %in% c("Negative", "Doublet") &
                                    dd$HTO_classification.comb2 == tag.var[i])] <- tag.var[i]
    }
    rm(i)
    
    
    var.hot <- paste("HTO_classification.r", kk, sep ="")
    dd[,var.hot] <- dd$HTO_classification
    rm(var.hot)
    
  }
  rm(kk)
  
  
  dd$classification.diff <- ifelse(dd$HTO_classification == dd$HTO_classification.tiger, 0, 1)
  
  dd$HTO_classification3 <- as.character(dd$HTO_classification)
  dd$HTO_classification3[which(dd$HTO_classification.tiger == "Doublet" & dd$classification.diff == 1)] <- "Doublet to Singlet"
  dd$HTO_classification3[which(dd$HTO_classification.tiger == "Negative" & dd$classification.diff == 1)] <- "Negative to Singlet"
  
  dd$HTO_classification3 <- factor(dd$HTO_classification3,
                                   levels = c(tag.var1, "Doublet", "Doublet to Singlet", "Negative to Singlet"))
  obj$HTO_classification3 <- dd$HTO_classification3
  
  end_time2 <- Sys.time()
  time_run2 <- end_time2 - start_time2
  time_run2 # Time difference of 10.02387 mins
  
  return(obj)
}

scDemultiplex_umap<-function(obj){
  tagnames<-rownames(obj)
  VariableFeatures(obj) <- tagnames
  obj <- ScaleData(obj)
  obj <- RunUMAP(obj, features=tagnames, slot="scale.data")
  return(obj)
}

scDemultiplex_plot<-function(obj, output_prefix){
  tagnames<-rownames(obj)

  g<-DimPlot(obj, reduction = "umap", group.by = "HTO_classification2")
  png(paste0(output_prefix, ".umap.class.png"), width=2000, height=1800, res=300)
  print(g)
  dev.off()
  
  width=max(1600, length(tagnames) * 1000)
  
  Idents(obj) <- "HTO_classification2"
  png(paste0(output_prefix, ".class.ridge.png"), width=width, height=max(1400, length(tagnames) * 300), res=300)
  print(RidgePlot(obj, assay = "HTO", features = tagnames, ncol = length(tagnames)))
  dev.off()
  
  png(paste0(output_prefix, ".class.dist.png"), width=width, height=max(1400, length(tagnames) * 500), res=300)
  rplot(obj, assay = "HTO", features = tagnames, identName="HTO_classification2")
  dev.off()
  
  if (length(tagnames) == 2) {
    png(paste0(output_prefix, ".class.point.png"), width=2000, height=1800, res=300)
    print(FeatureScatter(object = obj, feature1 = tagnames[1], feature2 = tagnames[2],group.by="HTO_classification2"))
    dev.off()
  }
  
  nwidth=ceiling(sqrt(length(tagnames)))
  nheight=ceiling(length(tagnames)/nwidth)
  png(paste0(output_prefix, ".umap.tag.png"), width=nwidth*1500, height=1500*nheight, res=300)
  g<-FeaturePlot(obj, features=tagnames, reduction = "umap")
  print(g)
  dev.off()
  
  hto_names <- c(tagnames, "Singlet", "Doublet", "Doublet to Singlet", "Negative to Singlet")
  cols=rep("gray", length(hto_names))
  names(cols)=hto_names
  cols[['Negative']]="blue"
  cols[["Doublet"]]="red"
  cols[["Doublet to Singlet"]]="green"
  cols[["Negative to Singlet"]]="gold"
  
  png(paste0(output_prefix, ".umap.all.png"), width=2000, height=1800, res=300)
  g<-DimPlot(obj, reduction = "umap", label=T, group.by="HTO_classification3", order=c("Negative", "Doublet", "Doublet to Singlet", "Negative to Singlet"))+
    scale_color_manual(values=cols)
  print(g)
  dev.off()
  
  return(obj)
}
