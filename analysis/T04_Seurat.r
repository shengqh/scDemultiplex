rm(list=ls()) 

root_dir="C:/projects/scratch/cqs/shengq2/papers/20221116_scrna_hto/"

#devtools::install_github("jabiru/tictoc")
library(Seurat)
library(tictoc)

sample="hto12"
for(sample in c("hto12", "pbmc")){
  sample_folder=paste0(root_dir, sample)
  setwd(sample_folder)
  
  rds_file=paste0(sample, ".obj.rds")
  obj=readRDS(rds_file)
  
  result_folder=paste0(sample_folder, "/HTODemux")
  if(!dir.exists(result_folder)){
    dir.create(result_folder)
  }
  setwd(result_folder)

  tic(paste0("starting ", sample, "...\n"))
  obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)
  toc1=toc()

  obj$HTODemux <- obj$HTO_classification
  obj$HTODemux[which(obj$HTO_classification.global == "Doublet")] <- "Doublet"
  
  saveRDS(list("HTODemux"=toc1), paste0(sample, ".HTODemux.tictoc.rds"))
  
  obj<-hto_plot(obj, paste0(sample, ".HTODemux"), group.by="HTODemux")

  saveRDS(obj, paste0(sample, ".HTODemux.rds"))
}
# ----------------------------------------------------------------

# 3. the heuristic classifier of MULTI-seq - R
# https://github.com/chris-mcginnis-ucsf/MULTI-seq
# https://satijalab.org/seurat/reference/multiseqdemux

sample="pbmc"
for(sample in c("hto12", "pbmc")){
  sample_folder=paste0(root_dir, sample)
  setwd(sample_folder)
  
  rds_file=paste0(sample, ".obj.rds")
  obj=readRDS(rds_file)
  
  result_folder=paste0(sample_folder, "/MULTIseqDemux")
  if(!dir.exists(result_folder)){
    dir.create(result_folder)
  }
  setwd(result_folder)
  
  tagnames<-rownames(obj)
  
  tic(paste0("starting ", sample, "...\n"))
  obj <- MULTIseqDemux(obj, assay = "HTO")
  toc1=toc()
  
  obj$MULTIseqDemux <- as.character(obj$MULTI_classification)
  obj$MULTIseqDemux[which(!obj$MULTI_classification %in% c(tagnames, "Negative"))] <- "Doublet"
  
  saveRDS(list("MULTIseqDemux"=toc1), paste0(sample, ".MULTIseqDemux.tictoc.rds"))
  
  obj<-hto_plot(obj, paste0(sample, ".MULTIseqDemux"), group.by="MULTIseqDemux")
  
  saveRDS(obj, paste0(sample, ".MULTIseqDemux.rds"))
}

