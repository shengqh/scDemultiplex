
rm(list=ls()) 

library(Seurat)
library(ggplot2)
library(patchwork)
library(hrbrthemes)

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

rdsfile = "C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/hto_samples_preparation/result/hto12.hto.rds" 
counts<-readRDS(rdsfile)
cutoff_startval<-1
output_prefix<-"hto12"

setwd("C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/cosine/")

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

groups<-unique(scaled$HTO_classification)
groups<-groups[order(groups)]
cl<-groups[2]
for(cl in unique(scaled$HTO_classification)){
  ingroup = rownames(scaled)[scaled$HTO_classification == cl]
  outgroup = rownames(scaled)[scaled$HTO_classification != cl]
  
  cosine.group<-data.frame(t(cosines[ingroup,]))
  cosine.median<-apply(cosine.group, 1, median)
  
  cos_df=data.frame("gcosine" = cosine.median, "group"=scaled$HTO_classification)
  
  class.sets=c(cl, "Negative", "Doublet")
  cos_df$group[!(cos_df$group %in% class.sets)]="Singlet"
  
  g1<-ggplot(cos_df, aes(x=group, y=gcosine)) + geom_violin() + 
    xlab("") + ylab("Cosine similarity") + 
    theme_bw()
  
  g2<-ggplot(cos_df, aes(x=gcosine)) + 
    geom_histogram(position = 'identity', bins=100) + 
    theme_bw()
  
  g3<-ggplot(cos_df, aes(x=gcosine, fill=group, color=group)) + 
    geom_histogram(alpha=0.5, position = 'identity', bins=100) + 
    theme_bw()
  
  g<-g1 + g2 + g3 + plot_layout(design = "AABBBB
AACCCC")  
  png(paste0(output_prefix, ".", cl, ".cosine.png"), width=3000, height=2000, res=300)
  print(g)
  dev.off()
}

