setwd("~/backup/lab/BT109VP")
rm(list=ls())
gc()
library(Seurat)
library(singlet)
library(dplyr)
int.seu <- qs::qread("int.seu.new.qs")
#######Normalization#####
gsc.seu <- SCTransform(gsc.seu,vars.to.regress =  c("percent.mt","S.Score","G2M.Score"),verbose = F)
gsc.seu <- RunPCA(gsc.seu,assay = "SCT")
ElbowPlot(gsc.seu)
gsc.seu <- RunHarmony(gsc.seu, group.by.vars = "orig.ident", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
gsc.seu <- RunUMAP(gsc.seu,reduction = "harmony",dims = 1:30)
gsc.seu <- RunTSNE(gsc.seu,dims = 1:30)
gsc.seu = FindNeighbors(gsc.seu,reduction = "harmony",dims=seq(30),verbose=F)

Sys.setenv(reticulate_python='/home/hzg/miniconda3/envs/r-reticulate/bin/python')
Sys.setenv(reticulate_python='/home/hzg/miniconda3/bin/python')
library(reticulate)
gsc.seu = FindClusters(gsc.seu,res=1,verbose=F,method = "igraph",algorithm = 1)	
DimPlot(gsc.seu,reduction = "umap",label = T)
all.markers <- FindAllMarkers(gsc.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25) 
top10.markers <- all.markers %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
###
gsc.seu <- qs::qread("gsc.seu.new.qs")
all.markers <- FindAllMarkers(gsc.seu, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.5) 
top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
qs::qsave(gsc.seu,"gsc.seu.new.qs")
###
gsc.seu <- qs::qread("~/backup/lab/BT109VP/gsc.seu.new.qs")
p1 <-DimPlot(gsc.seu,reduction = "umap",label = T)+ggplot2::ggtitle("Normalized")+DimPlot(gsc.seu,reduction = "umap",label = T,group.by = "Phase")+DimPlot(gsc.seu,reduction = "umap",label = T,group.by = "orig.ident")+ggplot2::ggtitle("Normalized")
########Unnormal Normal
DefaultAssay(gsc.seu) <-"RNA"
gsc.seu <- SCTransform(gsc.seu,vars.to.regress =  c("percent.mt"),verbose = F)
gsc.seu <- RunPCA(gsc.seu,assay = "SCT")
DefaultAssay(gsc.seu) <-"SCT"
gsc.seu <- RunUMAP(gsc.seu,reduction = "pca",dims = 1:30)
gsc.seu <- RunTSNE(gsc.seu,reduction = "pca",dims = 1:30)
gsc.seu = FindNeighbors(gsc.seu,reduction = "pca",dims=seq(30),verbose=F)
gsc.seu = FindClusters(gsc.seu,res=0.5,verbose=F,method = "igraph",algorithm = 1)	
p2<-DimPlot(gsc.seu,reduction = "umap",label = T)+ggplot2::ggtitle("Unnormalized")+DimPlot(gsc.seu,reduction = "umap",label = T,group.by = "Phase")+DimPlot(gsc.seu,reduction = "umap",label = T,group.by = "orig.ident")+ggplot2::ggtitle("Unnormalized")

#########SCT
DefaultAssay(gsc.seu) <-"RNA"
gsc.seu <- gsc.seu %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst",nfLeatures= 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50,verbose = F) %>% 
  RunUMAP(dims = 1:30)
gsc.seu = FindClusters(gsc.seu,res=0.6,verbose=F,method = "igraph",algorithm = 1)	
p3=DimPlot(gsc.seu,reduction = "umap",label = T)+ggplot2::ggtitle("VST")+DimPlot(gsc.seu,reduction = "umap",label = T,group.by = "Phase")+DimPlot(gsc.seu,reduction = "umap",label = T,group.by = "orig.ident")+ggplot2::ggtitle("Unnormalized")
p1/p2/p3
gsc.seu$cell_type <- Idents(gsc.seu)
qs::qsave(gsc.seu,file = "aa.qs")

aa <- FindMarkers(gsc.seu,ident.1 = "4",ident.2 = "10")
all.markers <- FindAllMarkers(gsc.seu, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.5) 
top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


int.seu <- qs::qread("gsc.seu.new.qs")

######
DefaultAssay(int.seu) <- "SCT"

int.seu$shank <- GetAssayData(int.seu)["SHANK3",]
int.seu$dlg <- GetAssayData(int.seu)["DLG4",]

int.seu$shank <- ifelse(int.seu$shank>0,"pos","neg")
int.seu$dlg <- ifelse(int.seu$dlg>0,"pos","neg")

meta <- int.seu@meta.data %>% group_by(orig.ident)
