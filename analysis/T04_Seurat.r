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
  
  result_folder=paste0(sample_folder, "/Seurat_HTODemux")
  if(!dir.exists(result_folder)){
    dir.create(result_folder)
  }
  setwd(result_folder)

  tic(paste0("starting ", sample, "...\n"))
  obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)
  toc1=toc()

  obj$Seurat_HTODemux <- obj$HTO_classification
  obj$Seurat_HTODemux[which(obj$HTO_classification.global == "Doublet")] <- "Doublet"
  
  saveRDS(list("HTODemux"=toc1), paste0(sample, ".Seurat_HTODemux.tictoc.rds"))
  
  obj<-hto_plot(obj, paste0(sample, ".Seurat_HTODemux"), group.by="Seurat_HTODemux")

  saveRDS(obj, paste0(sample, ".Seurat_HTODemux.rds"))
}
# ----------------------------------------------------------------

# 3. the heuristic classifier of MULTI-seq - R
# https://github.com/chris-mcginnis-ucsf/MULTI-seq
# https://satijalab.org/seurat/reference/multiseqdemux

sample="hto12"
for(sample in c("hto12", "pbmc")){
  sample_folder=paste0(root_dir, sample)
  setwd(sample_folder)
  
  rds_file=paste0(sample, ".obj.rds")
  obj=readRDS(rds_file)
  
  result_folder=paste0(sample_folder, "/Seurat_MULTIseqDemux")
  if(!dir.exists(result_folder)){
    dir.create(result_folder)
  }
  setwd(result_folder)
  
  tagnames<-rownames(obj)
  
  tic(paste0("starting ", sample, "...\n"))
  obj <- MULTIseqDemux(obj, assay = "HTO")
  toc1=toc()
  
  obj$Seurat_MULTIseqDemux <- as.character(obj$MULTI_classification)
  obj$Seurat_MULTIseqDemux[which(!obj$MULTI_classification %in% c(tagnames, "Negative"))] <- "Doublet"
  
  saveRDS(list("MULTIseqDemux"=toc1), paste0(sample, ".Seurat_MULTIseqDemux.tictoc.rds"))
  
  obj<-hto_plot(obj, paste0(sample, ".Seurat_MULTIseqDemux"), group.by="Seurat_MULTIseqDemux")
  
  saveRDS(obj, paste0(sample, ".Seurat_MULTIseqDemux.rds"))
}

