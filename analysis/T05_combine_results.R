rm(list=ls()) 

root_dir="C:/projects/scratch/cqs/shengq2/papers/20221116_scrna_hto/"

sample="pbmc"
for(sample in c("hto12", "pbmc")) {
  setwd(paste0(root_dir, sample))
  
  obj1<-readRDS(paste0("scDemultiplex/", sample, ".scDemultiplex.rds"))
  obj1.meta<-obj1@meta.data
  
  obj2.meta<-readRDS(paste0("HTODemux/", sample, ".HTODemux.rds"))@meta.data
  
  obj3.meta<-readRDS(paste0("MULTIseqDemux/", sample, ".MULTIseqDemux.rds"))@meta.data
  
  dd3 <- read.csv("GMM-demux/GMM_full.csv")
  dd3c <- read.csv("GMM-demux/GMM_full.config")
  gmm <- merge(dd3, dd3c, by.x = "Cluster_id", by.y = "X0", all.x = T, all.y = F)
  gmm$GMM_demux <- gmm$negative
  gmm$GMM_demux <- gsub(" ", "", gmm$GMM_demux)
  gmm$GMM_demux[which(is.na(gmm$negative))] <- "Negative"
  gmm$GMM_demux <- gsub("_", "-", gmm$GMM_demux)
  gmm$GMM_demux[which(! gmm$GMM_demux %in% c("Negative", rownames(obj1)))] <- "Doublet"
  rownames(gmm)<-gmm$X
  
  obj2.meta<-obj2.meta[rownames(obj1.meta),]
  obj3.meta<-obj3.meta[rownames(obj1.meta),]
  gmm<-gmm[rownames(obj1.meta),]
  
  obj1.meta$HTODemux = obj2.meta$HTODemux
  obj1.meta$MULTIseqDemux = obj3.meta$MULTIseqDemux
  obj1.meta$GMM_demux = gmm$GMM_demux
  
  write.csv(obj1.meta, file = paste0(sample, ".results.csv"), row.names = T)
  
  obj1@meta.data<-obj1.meta
  saveRDS(obj1, paste0(sample, ".results_obj.rds"))
}