##########Alghamri2021#########
rm(list = ls())
library(Seurat)
library(dplyr)
setwd("~/rna/sc/LGG/Alghamri2021/")
##Load the dataset you want to merge
path <- grep("GSM",dir(),value = T)
sample_names <- path
# Load and QC
scRNAlist <- list()
for (i in 1:6) {
  testdf <- Read10X(paste0("./",path[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = testdf$`Gene Expression`,project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id= sample_names[i])
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
}
for (i in 7:10) {
  testdf <- Read10X(paste0("./",path[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = testdf,project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id= sample_names[i])
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
}

# Simple Merge
Alghamri2021.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)

Alghamri2021.seu <- JoinLayers(Alghamri2021.seu)
table(Alghamri2021.seu$orig.ident)
Alghamri2021.seu[["percent.mt"]] <- PercentageFeatureSet(Alghamri2021.seu, pattern = "^MT-")
VlnPlot(Alghamri2021.seu, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
Alghamri2021.seu <- subset(Alghamri2021.seu, subset = percent.mt<25&nFeature_RNA >=250&nCount_RNA>=1000)
Alghamri2021.seu$patient <- Alghamri2021.seu$orig.ident
Alghamri2021.seu$patient <- gsub("GSM4610546|GSM4610547|GSM4610550|GSM4610551","7-16",Alghamri2021.seu$patient)
Alghamri2021.seu$patient <- gsub("GSM4610542|GSM4610543","3-20",Alghamri2021.seu$patient)
Alghamri2021.seu$patient <- gsub("GSM4610544|GSM4610545","4-30",Alghamri2021.seu$patient)
Alghamri2021.seu$patient <- gsub("GSM4610540|GSM4610541","3-18",Alghamri2021.seu$patient)
library(harmony)
standard10X = function(int.seu,nPCs=30,res=0.5,verbose=FALSE){
  int.seu = NormalizeData(int.seu)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu,group.by.vars="patient")
  int.seu = RunUMAP(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
Alghamri2021.seu <- standard10X(Alghamri2021.seu,verbose = T)
DimPlot(Alghamri2021.seu,group.by = "orig.ident",label = T)+DimPlot(Alghamri2021.seu,label = T)
aa <- FindAllMarkers(Alghamri2021.seu,logfc.threshold = 2,only.pos = T)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

qs::qsave(Alghamri2021.seu,file = "Alghamri2021-10x.seu.qs")
