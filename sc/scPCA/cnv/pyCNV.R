rm(list = ls())
gc()
library(Seurat)
int.seu <- qs::qread("int.seu.qs")
int.seu <- DietSeurat(int.seu)
# 第三文件基因注释文件
library(AnnoProbe)
DefaultAssay(int.seu) <- "RNA"
int.seu@assays$SCT <- NULL
geneInfor <- annoGene(rownames(int.seu), "SYMBOL", "human")
geneInfor <- geneInfor[with(geneInfor, order(chr, start)), c(1, 4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[, 1]), ]
geneInfor <- geneInfor[!geneInfor$chr %in% c("chrY", "chrX"), ] # nolint: object_name_linter. # nolint
int.seu <- subset(int.seu,features = geneInfor$SYMBOL)

library(SeuratDisk)
metadata <- int.seu@meta.data
write.csv(int.seu@meta.data, file = "./cnv/metadata.csv")
write.csv(geneInfor, file = "./cnv/gene.csv",row.names = F)

SaveH5Seurat(int.seu, filename = "./cnv/int.seu.h5Seurat")
Convert("./cnv/int.seu.h5Seurat", dest = "h5ad",overwrite = T)

####################
#####DownStream#####
####################
library(MetBrewer)
library(patchwork)
library(ggplot2)
library(dplyr)

seurat_colors <- as.character(sample(met.brewer("Klimt", 14),14))
colors_continuous <- as.character(rev(met.brewer("Hiroshige", 100)))
colors_scale <- as.character(rev(met.brewer("Johnson",100)))
sub_cols <- c("#D2C564","#CF5A79")

cnv.umap <- read.csv("cnv_umap.csv",header = T)[,-1]
cnv.score <- read.csv("./cnv/cnv_score.csv",header = T)
cnv.leiden <- read.csv("./cnv/cnv_leiden.csv",header = T)
int.seu <- qs::qread("./int.seu.qs")
int.seu$cnv.score <- cnv.score$cnv_score
int.seu$cnv.leiden <- as.factor(cnv.leiden$cnv_leiden)

VlnPlot(int.seu,features = "cnv.score",group.by = "cnv.leiden")
DimPlot(int.seu,group.by = "cnv.leiden",label = T)+
  FeaturePlot(int.seu,features = "cnv.score")

int.seu$tumor <- ifelse(int.seu$cnv.score>0.03,"tumor","normal")
DimPlot(int.seu,group.by = "tumor")
qs::qsave(int.seu,file = "scPCA_HGG.qs")

gsc.seu <- subset(int.seu,tumor=="tumor")

gsc.seu <- SCTransform(gsc.seu,vars.to.regress = "percent.mt")
gsc.seu = RunPCA(gsc.seu,verbose=F)
# gsc.seu = harmony::RunHarmony(gsc.seu,group.by.vars="orig.ident",max_iter=20)
gsc.seu <- RunUMAP(gsc.seu,reduction = "pca",dims = 1:30)
gsc.seu = FindNeighbors(gsc.seu,reduction = "pca",dims=seq(30),verbose=F)
gsc.seu = FindClusters(gsc.seu,res=0.5,verbose=F)
DimPlot(gsc.seu,reduction = "umap",label = T)
