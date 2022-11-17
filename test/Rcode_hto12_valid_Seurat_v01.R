
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


rdsfile <- "D:/consulting/043_SingleCell_Tiger/valid_data_pbmc_hto12/data/hto12.hto.rds"
output_prefix <- "hto12.HTO"

obj <- read_hto(rdsfile, output_prefix)

# An object of class Seurat 
# 12 features across 8193 samples within 1 assay 
# Active assay: HTO (12 features, 0 variable features)

tagnames <- rownames(obj[["HTO"]])

# [1] "HEK-A"  "HEK-B"  "HEK-C"  "THP1-A" "THP1-B" "THP1-C" "K562-A" "K562-B" "K562-C"
# [10] "KG1-A"  "KG1-B"  "KG1-C" 

data <- FetchData(object=obj, vars=tagnames)

# dim(data) # 8193   12

# head(data, 2)
#                      HEK-A    HEK-B     HEK-C    THP1-A    THP1-B    THP1-C
# CAGATCAAGTAGGCCA 0.4181826 3.868401 1.2540169 0.9765397 1.0019763 0.5616886
# CCTTTCTGTCGGATCC 0.0000000 0.000000 0.3058438 1.1303219 0.8666497 0.4480285
#                    K562-A    K562-B    K562-C    KG1-A     KG1-B     KG1-C
# CAGATCAAGTAGGCCA 4.066058 0.7629875 0.3361489 1.229882 1.1126684 0.9274255
# CCTTTCTGTCGGATCC 0.000000 4.8067667 0.7878399 1.104380 0.5393333 1.0002611

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
# 1257      260     6676

# table(obj2$HTO_classification)

obj2$HTO_classification2 <- obj2$HTO_classification
obj2$HTO_classification2[which(obj2$HTO_classification.global == "Doublet")] <- "Doublet"

# table(obj2$HTO_classification2)

# Doublet    HEK-A    HEK-B    HEK-C   K562-A   K562-B   K562-C    KG1-A 
#    1257      588      630      553      523      607      587      587 
#  KG1-B    KG1-C Negative   THP1-A   THP1-B   THP1-C 
#    509      582      260      457      502      551

data$HTO.classification.Seurat.HTODemux <- obj2$HTO_classification2


Idents(obj2) <- "HTO_classification2"
VariableFeatures(obj2) <- tagnames
obj2 <- ScaleData(obj2)
obj2 <- RunUMAP(obj2, features=tagnames, slot="scale.data")
DimPlot(obj2, reduction = "umap", group.by = "HTO_classification2")

png(paste0("../v01_Seurat/Seurat.HTODemux.hto12", ".umap.all.png"), width=2000, height=1800, res=300)
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

# Doublet    HEK-A    HEK-B    HEK-C   K562-A   K562-B   K562-C    KG1-A 
# 255      609      676      551      523      596      535      559 
# KG1-B    KG1-C Negative   THP1-A   THP1-B   THP1-C 
# 119      123     1981      518      541      607 


data$HTO.classification.Seurat.MULTIseqDemux <- obj3$MULTI_classification2

if(F){
  
  dd <- subset(data, select = c("HTO.classification.Seurat.HTODemux", 
                                "HTO.classification.Seurat.MULTIseqDemux"))
  
  write.csv(dd, file = "../../hto12_HTO_Seurat_results.csv")
  
  
}


Idents(obj3) <- "MULTI_classification2"
VariableFeatures(obj3) <- tagnames
obj3 <- ScaleData(obj3)
obj3 <- RunUMAP(obj3, features=tagnames, slot="scale.data")
DimPlot(obj3, reduction = "umap", group.by = "MULTI_classification2")

png(paste0("../v01_Seurat/Seurat.MULTIseqDemux.hto12", ".umap.all.png"), width=2000, height=1800, res=300)
g <- DimPlot(obj3, reduction = "umap", group.by = "MULTI_classification2")
print(g)
dev.off()

# ----------------------------------------------------------------





