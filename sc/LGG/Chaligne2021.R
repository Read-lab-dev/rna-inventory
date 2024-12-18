##########Chaligne2021#########
rm(list = ls())
library(Seurat)
library(dplyr)
path <- grep(".gz",dir("Chaligne2021/"),value = T)
tmp=as.data.frame(data.table::fread(paste0("./Chaligne2021/",path[1])))
rownames(tmp)<- tmp$GENE
tmp <- tmp[,-1]
exp <- tmp
for (i in 2:7) {
  tmp=as.data.frame(data.table::fread(paste0("./Chaligne2021/",path[i])))
  rownames(tmp)<- tmp$GENE
  tmp <- tmp[,-1]
  tmp <- as.data.frame(tmp)
  exp <- cbind(exp,tmp)
}
colSums((2^exp[,1:5]-1)*10) ###Check if it is TPM
colSums(exp[,1:5])

Chaligne2021.seu <- CreateSeuratObject(counts = exp)
Chaligne2021.seu$orig.ident <- stringr::str_split(colnames(Chaligne2021.seu),'[.]',simplify = T)[,1]
Chaligne2021.seu$orig.ident <- stringr::str_split(Chaligne2021.seu$orig.ident,'[_]',simplify = T)[,1]
table(Chaligne2021.seu$orig.ident)
Chaligne2021.seu[["percent.mt"]] <- PercentageFeatureSet(Chaligne2021.seu, pattern = "^MT")
VlnPlot(Chaligne2021.seu, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
library(harmony)
standard_SS2 = function(int.seu,nPCs=30,res=0.8,verbose=FALSE){
  LayerData(int.seu,layer = "data") = log2((LayerData(int.seu,layer = "counts")/10)+1)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu,group.by.vars="orig.ident")
  int.seu = RunUMAP(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
Chaligne2021.seu <- standard_SS2(Chaligne2021.seu,verbose = T)
DimPlot(Chaligne2021.seu,group.by = "orig.ident",label = T)+DimPlot(Chaligne2021.seu,label = T)
aa <- FindAllMarkers(Chaligne2021.seu,logfc.threshold = 2,only.pos = T)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Chaligne2021.seu <- subset(Chaligne2021.seu,seurat_clusters%in%c("5","9","8"),invert=T)

qs::qsave(Chaligne2021.seu,file = "Chaligne2021.seu.qs")
Chaligne2021.seu <- qs::qread("Chaligne2021.seu.qs")
