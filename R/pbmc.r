
rm(list=ls()) 

library(Seurat)
library(ggplot2)

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

source("C:/Users/sheng/Programs/scDemultiplex/R/utils.R")
source("C:/Users/sheng/Programs/scDemultiplex/R/scDemultiplex.R")

rdsfile = "C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/hto_samples_preparation/result/pbmc.hto.rds" 
counts<-readRDS(rdsfile)
cutoff_startval<-2
output_prefix<-"pbmc.HTO"

obj<-scDemultiplex_by_cutoff(counts, output_prefix, cutoff_startval)
obj<-scDemultiplex_umap(obj)

p.cut=0.0001
for (p.cut in c(0.0001, 0.00001, 0.000001)){
  obj2<-scDemultiplex_by_refine(obj, p.cut)
  obj2<-scDemultiplex_plot(obj2, paste0(output_prefix, ".p", p.cut))
}

