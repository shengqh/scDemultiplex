rm(list=ls()) 

root_dir="C:/projects/scratch/cqs/shengq2/papers/20221116_scrna_hto/"

#devtools::install_github("jabiru/tictoc")
library(scDemultiplex)
library(tictoc)

p.cut=0.001

sample="pbmc"
for(sample in c("hto12", "pbmc")){
  sample_folder=paste0(root_dir, sample)
  setwd(sample_folder)
  
  rds_file=paste0(sample, ".obj.rds")
  obj=readRDS(rds_file)

  result_folder=paste0(sample_folder, "/scDemultiplex")
  if(!dir.exists(result_folder)){
    dir.create(result_folder)
  }
  setwd(result_folder)

  final_file=paste0(sample, ".scDemultiplex.rds")
  
  output_prefix<-paste0(sample, ".HTO")

  tic(paste0("starting ", sample, "...\n"))
  obj<-demulti_cutoff(obj, output_prefix, 0)
  toc1=toc()
  obj<-demulti_refine(obj, p.cut)
  toc2=toc()

  saveRDS(list("cutoff"=toc1, "refine"=toc2), paste0(sample, ".scDemultiplex.tictoc.rds"))

  obj<-hto_plot(obj, paste0(output_prefix, ".cutoff"), group.by="scDemultiplex_cutoff")
  obj<-hto_plot(obj, paste0(output_prefix, ".refine_p", p.cut), group.by="scDemultiplex")

  saveRDS(obj, final_file)
}
