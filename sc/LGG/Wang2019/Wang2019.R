##########Wang2019#########
rm(list = ls())
library(Seurat)
library(dplyr)
setwd("~/rna/sc/LGG/Wang2019/")
##Load the dataset you want to merge
path <- grep("SF",dir(),value = T)
sample_names <- path
# Load and QC
scRNAlist <- list()
for (i in 1:length(path)) {
  testdf <- Read10X(paste0("./",path[i]),gene.column=1)
  scRNAlist[[i]] <- CreateSeuratObject(counts = testdf,project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id= sample_names[i])
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
}
# Simple Merge
Wang2019.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)

Wang2019.seu <- JoinLayers(Wang2019.seu)
table(Wang2019.seu$orig.ident)

Wang2019.seu[["percent.mt"]] <- PercentageFeatureSet(Wang2019.seu, pattern = "^MT-")
VlnPlot(Wang2019.seu, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
Wang2019.seu <- subset(Wang2019.seu, subset = percent.mt<20&nFeature_RNA >=250&nCount_RNA>=1000)
library(harmony)
standard10X = function(int.seu,nPCs=30,res=0.5,verbose=FALSE){
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
Wang2019.seu <- standard10X(Wang2019.seu,verbose = T)
DimPlot(Wang2019.seu,group.by = "orig.ident",label = T)+DimPlot(Wang2019.seu,label = T)
aa <- FindAllMarkers(Wang2019.seu,logfc.threshold = 2,only.pos = T)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

qs::qsave(Wang2019.seu,file = "Wang2019-10x.seu.qs")
