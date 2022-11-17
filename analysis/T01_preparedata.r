root_dir="C:/projects/scratch/cqs/shengq2/papers/20221116_scrna_hto/"

library(scDemultiplex)
library(zoo)
library("R.utils")
library(reshape2)
library(Matrix)

save_to_matrix<-function(counts, target_folder) {
  if(!dir.exists(target_folder)){
    dir.create(target_folder)
  }
  
  bar_file=paste0(target_folder, "/barcodes.tsv")
  writeLines(colnames(counts), bar_file)
  gzip(bar_file, overwrite=T)
  
  feature_file=paste0(target_folder, "/features.tsv")
  writeLines(rownames(counts), feature_file)
  gzip(feature_file, overwrite=T)
  
  matrix_file=paste0(target_folder, "/matrix.mtx")
  writeMM(counts, matrix_file)
  gzip(matrix_file, overwrite=T)
}

ignored_tags<-c("no_match","ambiguous","total_reads","bad_struct")
sample="pbmc"
for(sample in c("hto12", "pbmc")){
  sample_folder=paste0(root_dir, sample)
  if(!dir.exists(sample_folder)){
    dir.create(sample_folder)
  }
  setwd(sample_folder)
  
  hto_file=paste0("C:/projects/data/cqs/seurat_data/", sample, "_hto_mtx.rds")
  exp_file=paste0("C:/projects/data/cqs/seurat_data/", sample, "_umi_mtx.rds")

  hto=readRDS(hto_file)
  if(sample == "hto12"){
    hto=t(hto)
  }
  exp=readRDS(exp_file)
  
  common_cells=colnames(hto)[colnames(hto) %in% colnames(exp)]
  exp=exp[,common_cells]
  
  hto.exp <- CreateSeuratObject(counts = exp, min.features = 200)
  cells.valid<-colnames(hto.exp)
  hto<-hto[!(rownames(hto) %in% ignored_tags), cells.valid]

  counts=as.sparse(hto)
  
  save_to_matrix(counts=counts, target_folder="data")

  rds_file=paste0(sample, ".counts.rds")
  saveRDS(counts, rds_file)
  
  obj <- scDemultiplex:::read_hto(rds_file)
  obj<-hto_umap(obj)
  saveRDS(obj, paste0(sample, ".obj.rds"))
}
