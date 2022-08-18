
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

rdsfile = "C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/hto_samples_preparation/result/hto12.hto.rds" 
counts<-readRDS(rdsfile)
cutoff_startval<-0
output_prefix<-"hto12.HTO"

obj<-scDemultiplex_by_cutoff(counts, output_prefix, cutoff_startval)
obj<-scDemultiplex_umap(obj)

p.cuts<-c(0.001, 0.0001, 0.00001, 0.000001)
p.cuts<-c(0.001)

p.cut=p.cuts[1]
for (p.cut in p.cuts){
  obj2<-scDemultiplex_by_refine(obj, p.cut)
  obj2<-scDemultiplex_plot(obj2, paste0(output_prefix, ".p", p.cut))
}

