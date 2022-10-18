
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

#sample="hto12"
#cutoff_startval<-2
sample="pbmc"
cutoff_startval<-0

rdsfile = paste0("C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/hto_samples_preparation/result/", sample, ".hto.rds")
counts<-readRDS(rdsfile)
output_prefix<-paste0(sample, ".HTO")

setwd(r"(C:\projects\scratch\cqs\shengq2\papers\20210703_scrna_hto\scDemultiplex)")

obj<-scDemultiplex_by_cutoff(counts, output_prefix, cutoff_startval)
obj<-scDemultiplex_umap(obj)
obj<-scDemultiplex_plot(obj, paste0(output_prefix, ".cutoff"), group.by="scDemultiplex_cutoff")

p.cut=0.001
obj<-scDemultiplex_by_refine(obj, p.cut)
obj<-scDemultiplex_plot(obj, paste0(output_prefix, ".refine_p", p.cut), group.by="scDemultiplex_refine")

saveRDS(obj, paste0(sample, ".refined.rds"))
