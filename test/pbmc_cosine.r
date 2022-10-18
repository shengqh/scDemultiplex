
rm(list=ls()) 

library(Seurat)
library(ggplot2)
library(patchwork)

#devtools::install_github("shengqh/cutoff")
#install.packages("bbmle")
library(choisycutoff)
library(zoo)
library(reshape2)
library(gridExtra)
library(ggExtra)

library(TailRank) # dbb # Beta-Binomial Distribution

library(edgeR)
library(dirmult)
library(MGLM) # ddirmn
library(lsa)

source("C:/Users/sheng/Programs/scDemultiplex/R/utils.R")
source("C:/Users/sheng/Programs/scDemultiplex/R/scDemultiplex.R")

setwd("C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/cosine/")

data_df=data.frame(rdsfile = c("C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/hto_samples_preparation/result/pbmc.hto.rds",
                               "C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/hto_samples_preparation/result/hto12.hto.rds"), 
                   cutoff_startval = c(2, 1),
                   output_prefix = c("pbmc", "hto12"))

for(idx in c(1:nrow(data_df))){
  rdsfile=data_df$rdsfile[idx]
  cutoff_startval=data_df$cutoff_startval[idx]
  output_prefix=data_df$output_prefix[idx]
  
  counts<-readRDS(rdsfile)
  
  rdsfile = paste0(output_prefix, ".obj.rds")
  if (!file.exists(rdsfile)){
    obj<-scDemultiplex_by_cutoff(counts, output_prefix, cutoff_startval)
    obj<-scDemultiplex_umap(obj)
    saveRDS(obj, rdsfile)
  }else{
    obj<-readRDS(rdsfile)
  }
  
  scaled = obj@assays$HTO@scale.data
  
  cosine_rds<-paste0(output_prefix, ".cosine.rds")
  if(!file.exists(cosine_rds)){
    library(lsa)
    cosines<-cosine(scaled)
    saveRDS(cosines, cosine_rds)
  }else{
    cosines<-readRDS(cosine_rds)
  }
  
  scaled=data.frame(t(scaled))
  scaled$HTO_classification = obj$HTO_classification
  
  all_cosdf<-data.frame("cell"=rownames(scaled))
  groups<-unique(scaled$HTO_classification)
  groups<-groups[order(groups)]
  cl<-groups[2]
  for(cl in unique(scaled$HTO_classification)){
    ingroup = rownames(scaled)[scaled$HTO_classification == cl]
    outgroup = rownames(scaled)[scaled$HTO_classification != cl]
    
    cosine.group<-data.frame(t(cosines[ingroup,]))
    cosine.median<-apply(cosine.group, 1, median)
    #   
    #   cos_df=data.frame("gcosine" = cosine.median, "group"=scaled$HTO_classification)
    # 
    #   class.sets=c(cl, "Negative", "Doublet")
    #   cos_df$group[!(cos_df$group %in% class.sets)]="Singlet"
    # 
    #   g1<-ggplot(cos_df, aes(x=group, y=gcosine)) + geom_violin() + 
    #     xlab("") + ylab("Cosine similarity") + 
    #     theme_bw()
    #   
    #   g2<-ggplot(cos_df, aes(x=gcosine)) + 
    #     geom_histogram(position = 'identity', bins=100) + 
    #     theme_bw()
    # 
    #   g3<-ggplot(cos_df, aes(x=gcosine, fill=group, color=group)) + 
    #     geom_histogram(alpha=0.5, position = 'identity', bins=100) + 
    #     theme_bw()
    # 
    #   g<-g1 + g2 + g3 + plot_layout(design = "AABBBB
    # AACCCC")  
    #   png(paste0(output_prefix, ".", cl, ".cosine.png"), width=3000, height=2000, res=300)
    #   print(g)
    #   dev.off()
    
    all_cosdf[,cl] = cosine.median
  }
  
  #all_cosdf<-all_cosdf[, order(colnames(all_cosdf))]
  write.csv(all_cosdf, paste0(output_prefix, ".cosine.csv"), row.names=F)
  
  all_cosdf<-all_cosdf[,!(colnames(all_cosdf) %in% c("Negative", "Doublet"))]
  rownames(all_cosdf)<-all_cosdf$cell
  all_cosdf<-all_cosdf[c(2:ncol(all_cosdf))]
  
  tag_tb<-data.frame(table(obj$HTO_classification))
  colnames(tag_tb)<-c("HTO_classification", "single_tag")
  
  g<-DimPlot(obj, group.by="HTO_classification") + ggtitle("Single tag cutoff")
  for(cosine_cutoff in c(0.6, 0.65, 0.7, 0.75, 0.8)){
    cats <- apply(all_cosdf, 1, function(x){
      xnames<-colnames(all_cosdf)[x >= cosine_cutoff]
      if (length(xnames) == 0){
        return("Negative")
      }else if(length(xnames) > 1){
        return("Doublet")
      }else{
        return(xnames[1])
      }
    })
    
    hto_key = paste0("HTO_cosine_", cosine_cutoff)
    obj[[hto_key]]<-cats
    tb<-data.frame(table(cats))
    colnames(tb)<-c("HTO_classification", paste0("Cosine cutoff ", cosine_cutoff))
    
    tag_tb = merge(tag_tb, tb, by="HTO_classification", all.x=T, all.y=F)

    g2<-DimPlot(obj, group.by=hto_key) + ggtitle(paste0("Cosine cutoff ", cosine_cutoff))
    g<-g+g2
  }  
  g<-g +plot_layout(ncol=3)
  png(paste0(output_prefix, ".cosine.cutoff.png"), width=6600, height=4000, res=300)
  print(g)
  dev.off()
  
  write.csv(tag_tb, paste0(output_prefix, ".cosine.cutoff.csv"), row.names=F)
  saveRDS(obj@meta.data, paste0(output_prefix, ".cosine.meta.rds"))
}