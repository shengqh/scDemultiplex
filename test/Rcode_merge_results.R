

setwd("D:/consulting/043_SingleCell_Tiger/valid_data_pbmc_hto12")


# pbmc


dd <- read.csv("pbmc_HTO_demultiplex_results.csv")
names(dd) <- c("X", "HTO.classification.demultiplex.cutoff", "HTO.classification.demultiplex.refinement")

dd2 <- read.csv("pbmc_HTO_Seurat_results.csv")

dd3 <- read.csv("pbmc_demux/GMM_full.csv")
dd3c <- read.csv("pbmc_demux/GMM_full.config")
dd3 <- merge(dd3, dd3c, by.x = "Cluster_id", by.y = "X0", all.x = T, all.y = F)
dd3$HTO.classification.GMM.Demux <- dd3$negative
dd3$HTO.classification.GMM.Demux <- gsub(" ", "", dd3$HTO.classification.GMM.Demux)
dd3$HTO.classification.GMM.Demux[which(is.na(dd3$negative))] <- "Negative"
dd3$HTO.classification.GMM.Demux[which(! dd3$HTO.classification.GMM.Demux %in% 
                     c("Negative", "HTO_A", "HTO_B", "HTO_C", 
                       "HTO_D", "HTO_E", "HTO_F", "HTO_G", "HTO_H"))] <- "Doublet"
dd3$HTO.classification.GMM.Demux <- gsub("_", "-", dd3$HTO.classification.GMM.Demux)
dd3 <- subset(dd3, select = c("X", "HTO.classification.GMM.Demux"))


dd0 <- merge(dd, dd2, by = "X")

dd0 <- merge(dd0, dd3, by = "X")

write.csv(dd0, file = "pbmc_results_v02.csv", row.names = F)

rm(dd); rm(dd0); rm(dd2); rm(dd3); rm(dd3c)

# -------------

# hto12

dd <- read.csv("hto12_HTO_demultiplex_results.csv")
names(dd) <- c("X", "HTO.classification.demultiplex.cutoff", "HTO.classification.demultiplex.refinement")

dd2 <- read.csv("hto12_HTO_Seurat_results.csv")

dd3 <- read.csv("hto12_demux/GMM_full.csv")
dd3c <- read.csv("hto12_demux/GMM_full.config")
dd3c$negative <- gsub(" ", "", dd3c$negative)

dd3 <- merge(dd3, dd3c, by.x = "Cluster_id", by.y = "X0", all.x = T, all.y = F)
dd3$HTO.classification.GMM.Demux <- dd3$negative

dd3$HTO.classification.GMM.Demux[which(is.na(dd3$negative))] <- "Negative"
dd3$HTO.classification.GMM.Demux[which(! dd3$HTO.classification.GMM.Demux %in% 
                                         c("Negative", dd3c$negative[1:12]))] <- "Doublet"
dd3$HTO.classification.GMM.Demux <- gsub("_", "-", dd3$HTO.classification.GMM.Demux)
dd3 <- subset(dd3, select = c("X", "HTO.classification.GMM.Demux"))


dd0 <- merge(dd, dd2, by = "X")

dd0 <- merge(dd0, dd3, by = "X")

write.csv(dd0, file = "hto12_results_v02.csv", row.names = F)


rm(dd); rm(dd0); rm(dd2); rm(dd3); rm(dd3c)


