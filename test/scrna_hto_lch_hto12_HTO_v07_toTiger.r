
rm(list=ls()) 
outFile='scrna_hto'
parSampleFile1='D:/consulting/043_SingleCell_Tiger/20210703_scrna_hto_Tiger_example/20210703_scrna_hto/hto_samples_cutoff_all/result/fileList1.txt'
parSampleFile2='D:/consulting/043_SingleCell_Tiger/20210703_scrna_hto_Tiger_example/20210703_scrna_hto/hto_samples_cutoff_all/result/fileList2.txt'
parSampleFile3='D:/consulting/043_SingleCell_Tiger/20210703_scrna_hto_Tiger_example/20210703_scrna_hto/hto_samples_cutoff_all/result/fileList3.txt'
parFile1=''
parFile2=''
parFile3=''

#setwd(r"(C:\projects\scratch\cqs\shengq2\papers\20210703_scrna_hto\hto_samples_cutoff_all\result)")
setwd("D:/consulting/043_SingleCell_Tiger/20210703_scrna_hto_Tiger_example/lch_run/v06_beta_binomial/results/v02")

source("D:/consulting/043_SingleCell_Tiger/20210703_scrna_hto_Tiger_example/lch_run/v03/split_samples_utils_lch.r")

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

source("D:/consulting/043_SingleCell_Tiger/20210703_scrna_hto_Tiger_example/lch_run/v03/hto_find_cutoff_funs_Tiger.r")

# ----


library(edgeR)
library(dirmult)
library(MGLM) # ddirmn



# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

files_lines=read.table(parSampleFile1, sep="\t") # 14 2
files_lines$V1 <- gsub("C:/projects/scratch/cqs/shengq2/papers/20210703_scrna_hto/hto_samples_preparation/result/", 
                       "D:/consulting/043_SingleCell_Tiger/20210703_scrna_hto_Tiger_example/20210703_scrna_hto/hto_samples_preparation/result/", 
                       files_lines$V1)
files = split(files_lines$V1, files_lines$V2)

# files_lines
#          V2
# 1     COVID
# 2      CTRL
# 3  HYW_4701
# 4  HYW_4847
# 5  HYW_4880
# 6  HYW_4881
# 7  HYW_5386
# 8  HYW_5742
# 9  HYW_5755
# 10     LTS4
# 11     Pep2
# 12      SEB
# 13    hto12
# 14     pbmc


# ---

cutoffs_lines = read.table(parSampleFile2, sep="\t")
cutoffs = split(cutoffs_lines$V1, cutoffs_lines$V2)

# cutoffs
# $hto12
# [1] 2

# ---

params_lines = read.table(parSampleFile3, sep="\t")
params = split(params_lines$V1, params_lines$V2)
params$hto_ignore_exists = ifelse(params$hto_ignore_exists=="0", FALSE, TRUE)

# params
# $hto_ignore_exists
# [1] TRUE

# ---

# SEB # idx=14

# hto12 # idx=3

idx=3 

  fname = names(files)[idx] # "hto12"
  output_prefix = paste0(fname, ".HTO") # "hto12.HTO"
  output_file = paste0(output_prefix, ".csv") # "hto12.HTO.csv"
  
  # if(file.exists(output_file) & params$hto_ignore_exists){
  #   next
  # }
  # 
  
  rdsfile = files[[idx]] 
  cat(fname, ":", rdsfile, " ...\n")

  obj = read_hto(rdsfile, output_prefix)
  # An object of class Seurat 
  # 12 features across 8533 samples within 1 assay 
  # Active assay: HTO (12 features, 0 variable features)
  
  cutoff_point = ifelse(fname %in% names(cutoffs), as.numeric(cutoffs[[fname]]), 0)

  tagnames = rownames(obj[["HTO"]])
  # "HEK-A"  "HEK-B"  "HEK-C"  "THP1-A" "THP1-B" "THP1-C" "K562-A" "K562-B"
  # [9] "K562-C" "KG1-A"  "KG1-B"  "KG1-C"
  
  data <- FetchData(object=obj, vars=tagnames)
  
  # FetchData() 
  # Retrieves data (feature expression, PCA scores, metrics, etc.) for a set of cells in a Seurat object
  
  # head(data, 3)
  #                     HEK-A     HEK-B     HEK-C    THP1-A    THP1-B    THP1-C
  # CAGATCAAGTAGGCCA 0.421468 3.8810433 1.2583431 0.9853003 1.0129832 0.5672014
  # CCTTTCTGTCGGATCC 0.000000 0.0000000 0.3074411 1.1398333 0.8767359 0.4526629
  # CATATGGCATGGAATA 0.000000 0.3560652 1.2583431 1.0262294 0.6296225 0.7631339
  #                    K562-A    K562-B    K562-C    KG1-A     KG1-B     KG1-C
  # CAGATCAAGTAGGCCA 4.076104 0.7696758 0.3388700 1.237046 1.1219129 0.9345649
  # CCTTTCTGTCGGATCC 0.000000 4.8191601 0.7930300 1.111150 0.5450838 1.0077274
  # CATATGGCATGGAATA 4.797841 3.4457734 0.9606688 1.083951 2.5932053 1.0423947
  
  
  tagname=tagnames[1]  
  for (tagname in tagnames) {
    values=data[,tagname]
    values=values[values>0] # remove count = 0
    cat(paste0("get cutoff of ", tagname, " ...\n"))
    cutoff=get_cutoff(values, paste0(output_prefix, "_", tagname), cutoff_point)
    data[,paste0(tagname,"_pos")] = ifelse(data[,tagname]>cutoff, tagname, "Negative")
  }
  
  
  cat(paste0("get classification ...\n"))
  class_names=paste0(tagnames, "_pos")
  data$HTO_classification=unlist(apply(data, 1, function(x){
    xx=unique(x[class_names])
    if (length(xx) > 1){
      xx = xx[xx != "Negative"]
      if (length(xx) > 1) {
        return("Doublet")
      }
    }
    return(xx)
  }))
  
  
  data$HTO_classification.global=unlist(apply(data, 1, function(x){
    xx=unique(x[class_names])
    if (length(xx) > 1){
      xx = xx[xx != "Negative"]
      if (length(xx) > 1) {
        return("Doublet")
      }else{
        return("Singlet")
      }
    }
    return("Negative")
  }))
  
  
  obj[["HTO_classification"]] = data$HTO_classification
  obj[["HTO_classification.global"]] = data$HTO_classification.global
  
  # -----
  
  # table(data$HTO_classification.global)
  #
  # Doublet Negative  Singlet 
  #     597      646     7290 
  
  # table(data$HTO_classification)
  #
  # Doublet    HEK-A    HEK-B    HEK-C   K562-A   K562-B   K562-C    KG1-A 
  #     597      636      690      612      580      636      595      640 
  # KG1-B    KG1-C Negative   THP1-A   THP1-B   THP1-C 
  #   558      650      646      518      565      610 

# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

# data set in count
  
dd <- obj[["HTO"]]@counts
dd <- t(as.matrix(dd)) # 8533   12
dd <- as.data.frame(dd)

tag.var <- names(dd)
tag.var1 <- c(tag.var, "Negative")

dd$total.count <- apply(dd, 1, sum)

dd$id <- row.names(dd)


# ----

start_time2 <- Sys.time()

dd$HTO_classification.tiger <- dd$HTO_classification.r0 <- dd$HTO_classification <- data$HTO_classification
dd$HTO_classification.global.tiger <- dd$HTO_classification.global.r0 <- dd$HTO_classification.global <- data$HTO_classification.global

for(kk in 1:10){
  
  print(kk)
  # start_time1 <- Sys.time()

  out.alpha.est <- matrix(NA, 2, length(tag.var))
  for(i in 1:length(tag.var)){
    
    dds <- subset(dd, HTO_classification == tag.var[i], select = tag.var) # HTO_classification
    tag.sum <- apply(dds[,tag.var], 2, sum)
    tag.sum <- c(tag.sum[i], sum(tag.sum[-i])) # make it to k=2
    p.tu <- goodTuringProportions(tag.sum)
    dds2 <- data.frame(x1=dds[,tag.var[i]],
                       x2=apply(dds[,tag.var[-i]], 1, sum))
    fit <- dirmult(dds2, init = p.tu, epsilon=10^(-4), trace = FALSE)
    theta.est <- fit$theta
    alpha.est <- p.tu*(1-theta.est)/theta.est
    out.alpha.est[, i] <- alpha.est
    
    rm(dds); rm(tag.sum); rm(p.tu); rm(fit); rm(theta.est)
  }
  rm(i)

  # end_time1 <- Sys.time()
  # time_run1 <- end_time1 - start_time1
  # time_run1 # 

  # ----

  # start_time2 <- Sys.time()

  NN <- dd$total.count
  nn.tag <- length(tag.var) 
  for(j in 1:nn.tag){
  
    alpha.est <- out.alpha.est[,j]
    uu <- alpha.est[1]
    vv <- alpha.est[2]
  
    var.new <- paste(tag.var[j], ".pbb", sep = "")
    dd[,var.new] <- NA
    dd[,var.new] <- sapply(1:nrow(dd), function(x) pbb(dd[x, tag.var[j]], NN[x], uu, vv))
    dd[,var.new] <- ifelse(dd[,var.new] < 0.5, dd[,var.new], 1 - dd[,var.new])
  
    var.run <- paste(tag.var[j], ".pbb.r", kk, sep = "")
    dd[,var.run] <- dd[,var.new] 
    rm(var.run)
    
    var.new1 <- paste(var.new, ".ps.pvalue.fdr", sep = "")
    dd[,var.new1] <- p.adjust(dd[,var.new], method = "BH")
    
    var.run <- paste(var.new1, ".r", kk, sep = "")
    dd[,var.run] <- dd[,var.new1]
    rm(var.run)
    
    rm(var.new); rm(var.new1)
  }
  rm(j)

  # end_time2 <- Sys.time()
  # time_run2 <- end_time2 - start_time2
  # time_run2 # 

  # ----

  dd$HTO_classification.comb2 <- NA
  dd$HTO_classification.list2 <- NA
  #p.cut <- 0.025
  #varp <- paste(tag.var, ".pbb", sep = "")
  p.cut <- 0.001
  varp <- paste(tag.var, ".pbb.ps.pvalue.fdr", sep = "")

  for(i in 1:nrow(dd)){
    out <- dd[i,varp]
    out2 <- tag.var[which(out > p.cut)]
    dd$HTO_classification.list2[i] <- length(out2)
    dd$HTO_classification.comb2[i] <- paste(out2, collapse = "; ")
    rm(out); rm(out2)
  }
  rm(i)

  dd$HTO_classification.comb2[which(dd$HTO_classification.list2 == 0)] <- "Negative"

  var.hot <- paste("HTO_classification.comb2.r", kk, sep ="")
  dd[,var.hot] <- dd$HTO_classification.comb2
  rm(var.hot)
  
  var.hot <- paste("HTO_classification.list2.r", kk, sep ="")
  dd[,var.hot] <- dd$HTO_classification.list2
  rm(var.hot)
  
  # ----

  for(i in 1:length(tag.var)){
    dd$HTO_classification[which(dd$HTO_classification.list2 == 1 & 
                                dd$HTO_classification %in% c("Negative", "Doublet") &
                                dd$HTO_classification.comb2 == tag.var[i])] <- tag.var[i]
  }
  rm(i)


  var.hot <- paste("HTO_classification.r", kk, sep ="")
  dd[,var.hot] <- dd$HTO_classification
  rm(var.hot)

}
rm(kk)


dd$classification.diff <- ifelse(dd$HTO_classification == dd$HTO_classification.tiger, 0, 1)

save.image(file = "D:/consulting/043_SingleCell_Tiger/20210703_scrna_hto_Tiger_example/lch_run/v06_beta_binomial/hto12_HTO_pbb_v07.RData")

end_time2 <- Sys.time()
time_run2 <- end_time2 - start_time2
time_run2 # Time difference of 10.02387 mins



# ----

obj$HTO_classification2 = factor(obj$HTO_classification, levels=c(tag.var1, "Doublet"))

VariableFeatures(obj) <- tagnames
obj <- ScaleData(obj)
obj <- RunUMAP(obj, features=tagnames, slot="scale.data")
DimPlot(obj, reduction = "umap", group.by = "HTO_classification2")

png(paste0(output_prefix, ".v07.umap.class.png"), width=2000, height=1800, res=300)
g<-DimPlot(obj, reduction = "umap", group.by="HTO_classification2")
print(g)
dev.off()

# ---

width=max(1600, length(tagnames) * 1000)

Idents(obj) <- "HTO_classification2"
png(paste0(output_prefix, ".v07.class.ridge.png"), width=width, height=max(1400, length(tagnames) * 300), res=300)
print(RidgePlot(obj, assay = "HTO", features = tagnames, ncol = length(tagnames)))
dev.off()

png(paste0(output_prefix, ".v07.class.dist.png"), width=width, height=max(1400, length(tagnames) * 500), res=300)
rplot(obj, assay = "HTO", features = tagnames, identName="HTO_classification2")
dev.off()

if (length(tagnames) == 2) {
  png(paste0(output_prefix, ".v07.class.point.png"), width=2000, height=1800, res=300)
  print(FeatureScatter(object = obj, feature1 = tagnames[1], feature2 = tagnames[2],group.by="HTO_classification2"))
  dev.off()
}

# 

nwidth=ceiling(sqrt(length(tagnames)))
nheight=ceiling(length(tagnames)/nwidth)
png(paste0(output_prefix, ".v07.umap.tag.png"), width=nwidth*1500, height=1500*nheight, res=300)
g<-FeaturePlot(obj, features=tagnames, reduction = "umap")
print(g)
dev.off()

# 


dd$HTO_classification3 <- as.character(dd$HTO_classification)
dd$HTO_classification3[which(dd$HTO_classification.tiger == "Doublet" & dd$classification.diff == 1)] <- "Doublet to Singlet"
dd$HTO_classification3[which(dd$HTO_classification.tiger == "Negative" & dd$classification.diff == 1)] <- "Negative to Singlet"

dd$HTO_classification3 <- factor(dd$HTO_classification3,
                                 levels = c(tag.var1, "Doublet", "Doublet to Singlet", "Negative to Singlet"))
obj$HTO_classification3 <- dd$HTO_classification3
hto_names <- c(tag.var1, "Doublet", "Doublet to Singlet", "Negative to Singlet")
cols=rep("gray", length(hto_names))
names(cols)=hto_names
cols[['Negative']]="blue"
cols[["Doublet"]]="red"
cols[["Doublet to Singlet"]]="green"
cols[["Negative to Singlet"]]="gold"

png(paste0(output_prefix, "2.v07.umap.all.png"), width=2000, height=1800, res=300)
g<-DimPlot(obj, reduction = "umap", label=T, group.by="HTO_classification3", order=c("Negative", "Doublet", "Doublet to Singlet", "Negative to Singlet"))+
  scale_color_manual(values=cols)
print(g)
dev.off()



#--------------










