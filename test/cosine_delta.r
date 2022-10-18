
rm(list=ls()) 

library(Seurat)
library(ggplot2)
library(patchwork)

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
library(lsa)

source("C:/Users/sheng/Programs/scDemultiplex/R/utils.R")
source("C:/Users/sheng/Programs/scDemultiplex/R/scDemultiplex.R")

setwd("C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/cosine/")

theme_bw2 <- function () { 
  theme_bw() %+replace% 
    theme(
      strip.background = element_rect(fill = NA, colour = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),			
      axis.line = element_line(colour = "black", size = 0.5)
    )
}


#output_prefix<-"hto12"
output_prefix<-"pbmc"

obj<-readRDS(paste0(output_prefix, ".obj.rds"))

all_cosdf<-read.csv(paste0(output_prefix, ".cosine.csv"), header = T, row.names=1, stringsAsFactors = F)
all_cosdf<-all_cosdf[,!(colnames(all_cosdf) %in% c("Negative", "Doublet"))]

delta<-unlist(apply(all_cosdf, 1, function(x){
  y=x[order(x, decreasing = T)]
  y[1] - y[2]
}))

obj$cosine_delta<-delta

png(paste0(output_prefix, ".cosine_delta.png"), width=1600, height=1500, res=300)
print(FeaturePlot(obj, features="cosine_delta"))
dev.off()

tags<-colnames(all_cosdf)

meta<-readRDS(paste0(output_prefix, ".cosine.meta.rds"))

assert(all(rownames(meta) == rownames(all_cosdf)))

cutoff<-NULL

tag=tags[1]
for(tag in tags){
  other_tags<-tags[tags != tag]
  delta<-apply(all_cosdf, 1, function(x){
    x[tag] - max(x[other_tags])
  })
  
  cosdf<-all_cosdf
  cosdf$Delta<-unlist(delta)
  
  cosdf$HTO_classification = gsub('-', '.', meta$HTO_classification)
  cosdf$group = cosdf$HTO_classification
  class.sets=c(tag, "Negative", "Doublet")
  cosdf$group[!(cosdf$group %in% class.sets)]="Singlet"
  
  tag_cells<-rownames(cosdf)[cosdf$group == tag]
  doublet_cells<-rownames(cosdf)[cosdf$group == "Doublet"]
  negative_cells<-rownames(cosdf)[cosdf$group == "Negative"]
  
  subcells<-rownames(cosdf)[cosdf$group %in% c(tag, "Doublet", "Negative")]
  
  obj$group<-cosdf$group
  obj[[tag]]<-cosdf[,tag]
  
  tag_delta=paste0(tag, "_delta")
  obj[[tag_delta]]<-cosdf$Delta
  
  is_tag=paste0("is_", tag)
  
  g1<-ggplot(cosdf, aes_string(x=tag, y="Delta", fill="group", color="group")) +
    geom_point() + facet_grid(group~.) + theme_bw2() + NoLegend()
  
  g2<-ggplot(cosdf, aes_string(x="HTO_classification", y=tag)) + geom_violin() +
    xlab("") + ylab("Cosine similarity") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  g3<-ggplot(cosdf, aes_string(x=tag, fill="group", color="group")) +
    geom_histogram(alpha=0.5, position = 'identity', bins=100) +
    theme_bw()
  
  g4<-ggplot(cosdf, aes_string(x="Delta", fill="group", color="group")) +
    geom_histogram(alpha=0.5, position = 'identity', bins=100) +
    theme_bw()

  g5<-FeaturePlot(obj, cells=subcells, features=tag, split.by="group") + plot_annotation(title='Cosine similarity', theme = theme(plot.title = element_text(size = 18, hjust=0.5)))
  
  g6<-FeaturePlot(obj, cells=subcells, features=tag_delta, split.by="group") + plot_annotation(title='Cosine delta', theme = theme(plot.title = element_text(size = 18, hjust=0.5)))

  obj[[is_tag]]<-cosdf[,tag] > 0.5 & cosdf$Delta > 0.3
  subobj<-subset(obj, cells=subcells)
  g7<-DimPlot(subset(obj, cells=subcells), cells=subcells, group.by=is_tag, split.by="group") + ggtitle('Cosine>0.5 & delta>0.3')
  
  g<-g1 + g2 + g3 + g4 + g5 + g6 + g7 + plot_layout(design = "AABBCCDD
AAEEEEEE
AAFFFFFF
AAGGGGGG")
  png(paste0(output_prefix, ".", tag, ".cosine.png"), width=5000, height=3500, res=300)
  print(g)
  dev.off()
  
  p1<-DimPlot(subobj, cells.highlight = doublet_cells)+ NoLegend() + ggtitle("Doublet")
  p2<-DimPlot(subobj, cells.highlight = negative_cells) + NoLegend() + ggtitle("Negative")
  p3<-DimPlot(subobj, cells.highlight = tag_cells) + NoLegend() + ggtitle(tag)
  
  for (cosine_cut in c(0.4,0.5,0.6)){
    p<-p1+p2+p3
    for(delta_cut in c(0.1,0.2,0.3)){
      obj[[is_tag]]<-cosdf[,tag] > cosine_cut & cosdf$Delta > delta_cut
      new_doublet_cells<-colnames(obj)[obj$group == "Doublet" & !(unlist(obj[[is_tag]]))]
      new_negative_cells<-colnames(obj)[obj$group == "Negative" & !(unlist(obj[[is_tag]]))]
      new_tag_cells<-colnames(obj)[unlist(obj[[is_tag]])]
    
      p4<-DimPlot(subobj, cells.highlight = new_doublet_cells) + NoLegend() + ggtitle(paste0("Cosine>", cosine_cut, ", delta>", delta_cut))
      p5<-DimPlot(subobj, cells.highlight = new_negative_cells) + NoLegend() + ggtitle("")
      p6<-DimPlot(subobj, cells.highlight = new_tag_cells) + NoLegend() + ggtitle("")
      
      p<-p + p4 + p5 + p6
    }
    p = p + plot_layout(design = "ABC
DEF
GHI
JKL")    
    png(paste0(output_prefix, ".", tag, ".cosine_", cosine_cut, ".png"), width=3000, height=4000, res=300)
    print(p)
    dev.off()
  }
}

