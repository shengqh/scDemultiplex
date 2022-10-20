
rm(list=ls()) 

library(scDemultiplex)

#source("C:/Users/sheng/Programs/scDemultiplex/R/utils.R")
#source("C:/Users/sheng/Programs/scDemultiplex/R/scDemultiplex.R")

#sample="hto12"
#cutoff_startval<-0
sample="pbmc"
cutoff_startval<-0

rdsfile = paste0("C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/hto_samples_preparation/result/", sample, ".hto.rds")
counts<-readRDS(rdsfile)
output_prefix<-paste0(sample, ".HTO")

setwd(r"(C:\projects\scratch\cqs\shengq2\papers\20210703_scrna_hto\scDemultiplex)")

obj<-demulti_cutoff(counts, output_prefix, cutoff_startval)
obj<-hto_umap(obj)
obj<-hto_plot(obj, paste0(output_prefix, ".cutoff"), group.by="scDemultiplex_cutoff")

p.cut=0.001
obj<-demulti_refine(obj, p.cut)
obj<-hto_plot(obj, paste0(output_prefix, ".refine_p", p.cut), group.by="scDemultiplex_refine")

saveRDS(obj, paste0(sample, ".refined.rds"))
