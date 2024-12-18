##########Abdelfattah2022#########
rm(list = ls())
library(Seurat)
library(dplyr)
library(reticulate)
library(anndata)
library(Seurat) 
library(Matrix)
setwd("~/rna/sc/LGG/Abdelfattah2022/")
path <- grep("GSM",dir(),value = T)
scRNAlist <- list()
for (i in 1:length(path)) {
  testdf <- Read10X(paste0("./",path[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = testdf,project = path[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id= path[i])
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
}
# Simple Merge
Abdelfattah2022.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)

Abdelfattah2022.seu <- JoinLayers(Abdelfattah2022.seu)
table(Abdelfattah2022.seu$orig.ident)
Abdelfattah2022.seu$orig.ident <- ifelse(Abdelfattah2022.seu$orig.ident=="GSM5518638","MDAG11","MDAG3")
Abdelfattah2022.seu[["percent.mt"]] <- PercentageFeatureSet(Abdelfattah2022.seu, pattern = "^MT-")
VlnPlot(Abdelfattah2022.seu, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
Abdelfattah2022.seu <- subset(Abdelfattah2022.seu, subset = percent.mt<25&nFeature_RNA >=250&nCount_RNA>=1000)
library(harmony)
standard10X = function(int.seu,nPCs=30,res=0.3,verbose=FALSE){
  int.seu = NormalizeData(int.seu)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu,group.by.vars="orig.ident")
  int.seu = RunUMAP(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
Abdelfattah2022.seu <- standard10X(Abdelfattah2022.seu,verbose = T)
DimPlot(Abdelfattah2022.seu,group.by = "orig.ident",label = T)+DimPlot(Abdelfattah2022.seu,label = T)
aa <- FindAllMarkers(Abdelfattah2022.seu,logfc.threshold = 2,only.pos = T)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

qs::qsave(Abdelfattah2022.seu,file = "Abdelfattah2022-10x.seu.qs")
