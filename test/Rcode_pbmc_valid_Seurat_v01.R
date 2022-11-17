
# valid data

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


setwd("D:/consulting/043_SingleCell_Tiger/valid_data_pbmc_hto12/results/v01")

source("D:/consulting/043_SingleCell_Tiger/20210703_scrna_hto_Tiger_example/lch_run/v03/split_samples_utils_lch.r")

source("D:/consulting/043_SingleCell_Tiger/20210703_scrna_hto_Tiger_example/lch_run/v03/hto_find_cutoff_funs_Tiger.r")

# ---------------------------------------------------------------------------------

rdsfile <- "D:/consulting/043_SingleCell_Tiger/valid_data_pbmc_hto12/data/pbmc.hto.rds"
output_prefix <- "pbmc.HTO"

obj <- read_hto(rdsfile, output_prefix)

# An object of class Seurat 
# 12 features across 9779 samples within 1 assay 
# Active assay: HTO (8 features, 0 variable features)

tagnames <- rownames(obj[["HTO"]])

# [1] "HTO-A" "HTO-B" "HTO-C" "HTO-D" "HTO-E" "HTO-F" "HTO-G" "HTO-H"

data <- FetchData(object=obj, vars=tagnames)

# dim(data) # 9779    8

# head(data, 2)
#                      HTO-A     HTO-B     HTO-C     HTO-D     HTO-E     HTO-F     HTO-G
# TAGTTGGGTCATACTG 0.3034395 0.4876152 0.2870507 0.6829806 0.6266393 0.1092183 0.2563741
# CTACACCGTACCGTAT 0.0000000 0.9225196 0.4873769 0.2599768 0.6881980 3.0793464 0.4306972
#                     HTO-H
# TAGTTGGGTCATACTG 3.341474
# CTACACCGTACCGTAT 0.573389

# ---------------------------------------------------------------------------------

# 1. the GMM-Demux classifier - Python
# 2. the heuristic classifier of Seurat - R
# 3. the heuristic classifier of MULTI-seq - R
# 4. the model-based classifier demuxEM - Python
# 5. a human-supervised classifier

# 1. the GMM-Demux classifier - Python
# https://github.com/CHPGenetics/GMM-Demux

# 4. the model-based classifier demuxEM
# https://github.com/chris-mcginnis-ucsf/MULTI-seq


# -------------------------------------------------------------

# 2. the heuristic classifier of Seurat
# https://satijalab.org/seurat/articles/hashing_vignette.html

library(Seurat)

obj2 <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)

# table(obj2$HTO_classification.global)

# Doublet Negative  Singlet 
#    2196       94     7489

# table(obj2$HTO_classification)

obj2$HTO_classification2 <- obj2$HTO_classification
obj2$HTO_classification2[which(obj2$HTO_classification.global == "Doublet")] <- "Doublet"

# table(obj2$HTO_classification2)

# Doublet    HTO-A    HTO-B    HTO-C    HTO-D    HTO-E    HTO-F    HTO-G    HTO-H 
#    2196      967     1223     1003      780      803      804      819     1090 
# Negative 
#      94  

data$HTO.classification.Seurat.HTODemux <- obj2$HTO_classification2


Idents(obj2) <- "HTO_classification2"
VariableFeatures(obj2) <- tagnames
obj2 <- ScaleData(obj2)
obj2 <- RunUMAP(obj2, features=tagnames, slot="scale.data")
DimPlot(obj2, reduction = "umap", group.by = "HTO_classification2")

png(paste0("../v01_Seurat/Seurat.HTODemux.pbmc", ".umap.all.png"), width=2000, height=1800, res=300)
g <- DimPlot(obj2, reduction = "umap", group.by = "HTO_classification2")
print(g)
dev.off()

# ----------------------------------------------------------------

# 3. the heuristic classifier of MULTI-seq - R
# https://github.com/chris-mcginnis-ucsf/MULTI-seq
# https://satijalab.org/seurat/reference/multiseqdemux

library(Seurat)

obj3 <- MULTIseqDemux(obj)

obj3$MULTI_classification2 <- as.character(obj3$MULTI_classification)
obj3$MULTI_classification2[which(!obj3$MULTI_classification2 %in% c(tagnames, "Negative"))] <- "Doublet"

# table(obj3$MULTI_classification2)
#
# Doublet    HTO-A    HTO-B    HTO-C    HTO-D    HTO-E    HTO-F    HTO-G    HTO-H 
#    1767      936     1253     1050      827      814      801      665     1089 
# Negative 
#     577

data$HTO.classification.Seurat.MULTIseqDemux <- obj3$MULTI_classification2

if(F){
  
  dd <- subset(data, select = c("HTO.classification.Seurat.HTODemux", 
                                "HTO.classification.Seurat.MULTIseqDemux"))
  
  write.csv(dd, file = "../../pbmc_HTO_Seurat_results.csv")
  
  
}


Idents(obj3) <- "MULTI_classification2"
VariableFeatures(obj3) <- tagnames
obj3 <- ScaleData(obj3)
obj3 <- RunUMAP(obj3, features=tagnames, slot="scale.data")
DimPlot(obj3, reduction = "umap", group.by = "MULTI_classification2")

png(paste0("../v01_Seurat/Seurat.MULTIseqDemux.pbmc", ".umap.all.png"), width=2000, height=1800, res=300)
g <- DimPlot(obj3, reduction = "umap", group.by = "MULTI_classification2")
print(g)
dev.off()

# ----------------------------------------------------------------





