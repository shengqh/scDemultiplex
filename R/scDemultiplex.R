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
library(doParallel)
library(tictoc)

check_mc_cores<-function(mc.cores) {  
  # if(.Platform$OS.type == "windows") {
  #   mc.cores=1
  # }else{
  #mc.cores=min(parallel::detectCores() - 1, max(1, mc.cores))
  # }
  mc.cores=min(parallel::detectCores(), max(1, mc.cores))
  return(mc.cores)
}

do_cutoff<-function(tagname, data, output_prefix=NULL, cutoff_startval=0){
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
  #don't use message in function called in thread, it might bring trouble.
  print(paste0("get cutoff of ", tagname, " ..."))
  if(is.null(output_prefix)){
    cur_prefix = NULL  
  }else{
    cur_prefix = paste0(output_prefix, ".", tagname)
  }

  ntags = ncol(data)

  cutoff=get_cutoff(tagname, values, ntags, cur_prefix, cur_cutoff)
  return(cutoff)
}

do_cutoff_parallel<-function(tagnames, data, output_prefix, cutoff_startval, mc.cores){
  if (is_windows() & (mc.cores > 1)) {
    print(paste("using", mc.cores, "threads in", .Platform$OS.type, " by parLapply."))
    cl <- makeCluster(mc.cores)  
    registerDoParallel(cl)  
    #https://stackoverflow.com/questions/12023403/using-parlapply-and-clusterexport-inside-a-function
    clusterExport(cl,list('zoo','rollapply', 'my_startval', 'my_cutoff', 'get_cutoff', "my_em", 'do_cutoff','data',"output_prefix","cutoff_startval"), envir=environment())
    system.time(
      results<-unlist(parLapply(cl,tagnames,fun=do_cutoff, data, output_prefix, cutoff_startval))
    )
    stopCluster(cl)
  }else{
    print(paste("using", mc.cores, "threads in", .Platform$OS.type, " by mclapply."))
    system.time(
      results<-unlist(parallel::mclapply(tagnames, do_cutoff, data = data, output_prefix = output_prefix, cutoff_startval = cutoff_startval, mc.cores=mc.cores))
    )
  }
  return(results)
}

##################################

#' @title demulti_cutoff:
#' @description Function that performs the initial cell classification (Singlet, Doublet, Negative).
#'  
#' @param counts = SEURAT OBJECT with the data to be analyzed.
#'
#' @param output_prefix = CHARACTER prefix to be put on the output file names (default NULL).
#' @param cutoff_startval = INTEGER Start value used to estimate cutoff (default 0).
#' @param mc.cores = INTEGER of how many cores to run in parallel (default 1).
#' @param cutoff_list = NAMED LIST cutoffs for each tagname (default NULL). If the cutoff_list is set, classification will be performed based on those predefined cutoffs.
#'
#' @returns Seurat object with initial characterizations
#' 
#' @examples 
#' #Load in and prepare the data
#' hto_counts = data.frame(fread("https://raw.githubusercontent.com/Oshlack/hashtag-demux-paper/main/data/batch1_c1_hto_counts.csv"), row.names=1, check.names=F)
#' #As an example, we sample 4000 cells to speed up
#' hto_counts <- hto_counts[, sample(colnames(hto_counts), size=4000, replace=F)]
#' #change to Seurat object
#' obj <- CreateSeuratObject(counts = hto_counts, assay="HTO")
#' obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
#' DefaultAssay(object = obj) <- "HTO"
#' tagnames=rownames(obj[["HTO"]])
#' ntags=length(tagnames)
#' 
#' #Normalize and run UMAP
#' obj<-ScaleData(obj, features=tagnames, verbose=FALSE)
#' obj<-RunPCA(obj, features=tagnames, approx=FALSE)
#' obj<-RunUMAP(obj, features=tagnames, slot="scale.data")
#' 
#' #Running demulti_cutoff
#' obj = demulti_cutoff(obj, mc.cores=ntags)
#' 
#' #Table of cells by Single HTO/Doublet/Negative
#' print(kable(table(obj$scDemultiplex_cutoff)))
#' 
#' #Table of Singlet/Foublet/Negative
#' print(kable(table(obj$scDemultiplex_cutoff.global)))
#' 
#' #Plotting the results of demulti_cutoff in UMAP format
#' g<-DimPlot(obj, reduction = "umap", group.by="scDemultiplex_cutoff")
#' print(g)
#' 
#' @export
demulti_cutoff<-function(counts, output_prefix=NULL, cutoff_startval=0, mc.cores=1, cutoff_list=NULL){
  if(is(counts,"Seurat")){
    obj=counts
  }else{
    obj<-get_object(counts)
  }
  
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
    if(!is.null(output_prefix)){
      for (tagname in tagnames) {
        cutoff = cutoff_list[tagname]
        values = data[,tagname]
        png(paste0(output_prefix, ".", tagname, ".cutoff.png"), width=2000, height=1600, res=300)
        hist(values,200,F,xlab="concentration",ylab="density", main=NULL,col="grey")
        lines(density(values),lwd=1.5,col="blue")
        abline(v=cutoff,lwd=1.5,col="brown")
        dev.off()
      }
    }
  }else{
    if(!require("choisycutoff")){
      BiocManager::install('shengqh/cutoff')
      library(choisycutoff)
    }
    mc.cores<-check_mc_cores(mc.cores)
    cutoff_list = do_cutoff_parallel(
      tagnames = tagnames, 
      data = data, 
      output_prefix = output_prefix, 
      cutoff_startval = cutoff_startval, 
      mc.cores=mc.cores )
    names(cutoff_list) = tagnames
  }
  if(!is.null(output_prefix)){
    saveRDS(cutoff_list, paste0(output_prefix, ".cutoff_list.rds"))
  }

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
  
##########################################

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

do_estimate_alpha_parallel<-function(tagnames, taglist, mc.cores){
  if (is_windows() & (mc.cores > 1)) {
    cat("using", mc.cores, "threads in", .Platform$OS.type, " by parLapply.\n")
    cl <- makeCluster(mc.cores)  
    registerDoParallel(cl)  
    clusterExport(cl,list('estimate_alpha', 'goodTuringProportions', 'dirmult', 'taglist'), envir=environment())
    system.time(
      results<-parLapply(cl,tagnames,fun=estimate_alpha, taglist)
    )
    stopCluster(cl)
  }else{
    cat("using", mc.cores, "threads in", .Platform$OS.type, " by mclapply.\n")
    system.time(
      results <- parallel::mclapply(tagnames, estimate_alpha, taglist=taglist, mc.cores=mc.cores)
    )
  }
  return(results)
}

#If at least 3 singlets moved from one tag (A) to another tag (B), we call it one cross assignment of B.
#If there are at least two cross assignment of B, the refinement stops.
should_stop<-function(begin_calls, refined_calls, min_singlet_cross_assigned=3, min_tag_cross_assigned=2){
  move_tb = table(begin_calls, refined_calls)
  move_tb_2 = move_tb[!(rownames(move_tb) %in% c("Doublet", "Negative")), !(colnames(move_tb) %in% c("Doublet", "Negative"))]
  
  move_from_other_singlets = unlist(apply(move_tb_2, 2, function(x){ sum(x >= min_singlet_cross_assigned) }))
  #since there is one value from self to self, using ">" instead of ">="
  result = any(move_from_other_singlets > min_tag_cross_assigned)
  return(result)
}

######################################

#' @title demulti_refine:
#' @description Function that refines the cell classification iterations performed in demulti_cutoff and calculates a pvalue representing the likelihood of a tag representing a doublet (two cells) or negative (zero cells).
#'
#' @param obj = SEURAT OBJECT with the data to be analyzed.
#'
#' @param output_prefix = CHARACTER prefix to be put on the output file names (default NULL).
#' @param p.cut = NUMERIC value between 0 and 1 to use as the pvalue cut off for significance (default 0.001).
#' @param iterations = INTEGER of how many cell classification iterations to go through (default 10).
#' @param init_column = CHARACTER name of the column with the initial classifications. If run directly after demulti_cutoff, leave as default (default "scDemultiplex_cutoff").
#' @param mc.cores = INTEGER of how many cores to run in parallel (default 1).
#' @param refine_negative_doublet_only = LOGICAL if want only negative cells reported, set as TRUE (default FALSE).
#' @param min_singlet_cross_assigned = INTEGER of how many cells are allowed to shift from singlet to doublet before halting the analysis (default 3).
#' @param min_tag_cross_assigned = INTEGER of how many cells can share tags before halting the analysis (default 2).
#'
#' @returns Seurat object with refined cell classifications
#' 
#' @examples
#' #This occurs after the example for demulti_cutoff
#' #object from the end of that example used as 'obj'
#' 
#' #Running demulti_refine
#' obj<-demulti_refine(obj, mc.cores=ntags)
#' 
#' #Table of Single HTO/Doublet/Negative
#' print(kable(table(obj$scDemultiplex)))
#' 
#' #Table of Singlet/Doublet/Negative
#' print(kable(table(obj$scDemultiplex.global)))
#' 
#' #Plotting the results of demulti_refine in a UMAP format
#' g<-DimPlot(obj, reduction = "umap", group.by="scDemultiplex")
#' print(g)
#' 
#' @export
demulti_refine<-function(obj, output_prefix=NULL, p.cut=0.001, iterations=10, init_column="scDemultiplex_cutoff", mc.cores=1, refine_negative_doublet_only=FALSE, min_singlet_cross_assigned=3, min_tag_cross_assigned=2){
  mc.cores<-check_mc_cores(mc.cores)

  dd <- obj[["HTO"]]@counts
  dd <- t(as.matrix(dd)) # 8193   12
  dd <- as.data.frame(dd)
  
  tag.var <- names(dd)

  init_classification = unlist(obj[[init_column]])
  classified_tags = unique(init_classification)

  if(!all(tag.var %in% classified_tags)){
    missed = tag.var[! tag.var %in% classified_tags]
    print(paste0("Warning: missing tags in the init classification, use init classification as final result: ", paste0(missed, collapse = ", ")))
    obj$scDemultiplex = init_classification
    obj$scDemultiplex.global = ifelse(obj$scDemultiplex %in% c("Negative", "Doublet"), as.character(obj$scDemultiplex), "Singlet")
  }else{
    nn.tag <- length(tag.var)
    dd$total.count <- apply(dd, 1, sum)
    NN <- dd$total.count
    
    start_time2 <- Sys.time()

    dd$HTO_classification <- init_classification
    df<-data.frame(table(dd$HTO_classification))
    df$Iteration=0

    hc<-dd[,"HTO_classification", drop=FALSE]
    colnames(hc)="0"

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
      out.alpha.est <- do_estimate_alpha_parallel(
        tagnames = names(taglist), 
        taglist=taglist, 
        mc.cores=mc.cores)
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
        print("  no change anymore, stop.")
        break
      }
      
      cur_df<-data.frame(table(dd$HTO_classification))
      cur_df$Iteration=kk   
      df<-rbind(df, cur_df) 
      
      cur_hc<-dd[,"HTO_classification", drop=FALSE]
      colnames(cur_hc)=kk
      hc<-cbind(hc, cur_hc) 

      if(kk > 1 & should_stop(lastClassification, dd$HTO_classification, min_singlet_cross_assigned, min_tag_cross_assigned)){
        dd$HTO_classification = lastClassification
        print("  too many singlets shifted from multiple tags to another same tag, stop.")
        break
      }
    }
    rm(kk)
    
    if(!is.null(output_prefix)){
      mdf<-dcast(df, Var1~Iteration, value.var="Freq")
      write.csv(mdf, paste0(output_prefix, ".iteration.csv"), row.names=FALSE)
      write.csv(hc, paste0(output_prefix, ".iteration.detail.csv"), row.names=TRUE)
    }
    
    end_time2 <- Sys.time()
    time_run2 <- end_time2 - start_time2
    time_run2
    
    # ----
    
    obj$scDemultiplex = factor(dd$HTO_classification, levels=c(tag.var, "Negative", "Doublet"))
    obj$scDemultiplex.global = ifelse(obj$scDemultiplex %in% c("Negative", "Doublet"), as.character(obj$scDemultiplex), "Singlet")
  }

  return(obj)
}

###########################################

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
