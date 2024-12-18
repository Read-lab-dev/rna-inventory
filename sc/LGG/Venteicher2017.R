##########Venteicher2017#########
rm(list = ls())
library(Seurat)
library(dplyr)
exp <- data.table::fread("./mis/GSE89567_IDH_A_processed_data.txt.gz")
exp <- as.data.frame(exp)
rownames(exp)<- gsub("'","",exp$V1)
exp <- exp[,-1]
colSums((2^exp[,1:5]-1)*10) ###Check if it is TPM

Venteicher2017.seu <- CreateSeuratObject(counts = exp)
Venteicher2017.seu$orig.ident <- stringr::str_split(colnames(Venteicher2017.seu),'-',simplify = T)[,1]
Venteicher2017.seu$orig.ident <- stringr::str_split(Venteicher2017.seu$orig.ident,'_',simplify = T)[,1]
table(Venteicher2017.seu$orig.ident)
Venteicher2017.seu$orig.ident[Venteicher2017.seu$orig.ident=="mgh103"]= "MGH103"
Venteicher2017.seu$orig.ident[Venteicher2017.seu$orig.ident%in%c("MGH107pos","MGH107neg")]= "MGH107"
Venteicher2017.seu$orig.ident[-grep("MGH",Venteicher2017.seu$orig.ident)] <- paste0("MGH",Venteicher2017.seu$orig.ident[-grep("MGH",Venteicher2017.seu$orig.ident)])

Venteicher2017.seu[["percent.mt"]] <- PercentageFeatureSet(Venteicher2017.seu, pattern = "^MT")
VlnPlot(Venteicher2017.seu, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
library(harmony)
standard_SS2 = function(int.seu,nPCs=30,res=0.8,verbose=FALSE){
  LayerData(int.seu,layer = "data") = LayerData(int.seu,layer = "counts")
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu,group.by.vars="orig.ident")
  int.seu = RunUMAP(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
Venteicher2017.seu <- standard_SS2(Venteicher2017.seu,verbose = T)
DimPlot(Venteicher2017.seu,group.by = "orig.ident",label = T)+DimPlot(Venteicher2017.seu,label = T)
aa <- FindAllMarkers(Venteicher2017.seu,logfc.threshold = 2,only.pos = T)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Venteicher2017.seu <- subset(Venteicher2017.seu,seurat_clusters%in%c("5","9","8"),invert=T)

qs::qsave(Venteicher2017.seu,file = "Venteicher2017.seu.qs")
Venteicher2017.seu <- qs::qread("Venteicher2017.seu.qs")
