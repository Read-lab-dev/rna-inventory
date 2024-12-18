#Integrate:smart-seq
library(Seurat)
library(harmony)
library(dplyr)
Tirosh2016.seu <- qs::qread("Tirosh2016.seu.qs")
colSums((2^LayerData(Tirosh2016.seu,layer = "counts")[,1:5]-1)*10) ###Check if it is TPM
Venteicher2017.seu <- qs::qread("Venteicher2017.seu.qs")
colSums((2^LayerData(Venteicher2017.seu,layer = "counts")[,1:5]-1)*10) ###Check if it is TPM
Spitzer2024.seu <- qs::qread("Spitzer2024.seu.qs")
colSums((2^LayerData(Venteicher2017.seu,layer = "counts")[,1:5]-1)*10) ###Check if it is TPM
Chaligne2021.seu <- qs::qread("Chaligne2021.seu.qs")
colSums((2^LayerData(Chaligne2021.seu,layer = "counts")[,1:5]-1)*10) ###Check if it is TPM
LayerData(Chaligne2021.seu,layer = "counts") <- log2((LayerData(Chaligne2021.seu,layer = "counts")/10)+1)

ss2.seu <- merge(Tirosh2016.seu,c(Venteicher2017.seu,Spitzer2024.seu,Chaligne2021.seu))
ss2.seu <- JoinLayers(ss2.seu)
table(ss2.seu$orig.ident)
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
ss2.seu <- standard_SS2(ss2.seu,verbose = T)
ss2.seu$platform <- "Smart-seq2"
DimPlot(ss2.seu,group.by = "orig.ident",label = T)+DimPlot(ss2.seu,label = T)

########10X###############
Spitzer2024.10x.seu <- qs::qread("Spitzer2024-10x.seu.qs") 
Wang2019.seu <- qs::qread("./Wang2019/Wang2019-10x.seu.qs") 
Alghamri2021.seu <- qs::qread("./Alghamri2021/Alghamri2021-10x.seu.qs")
Johnson2021.seu <- qs::qread("./Johnson2021/Johnson2021-10x.seu.qs")
Abdelfattah2022.seu <- qs::qread("./Abdelfattah2022/Abdelfattah2022-10x.seu.qs")
x10.seu <- merge(Spitzer2024.10x.seu,c(Wang2019.seu,Alghamri2021.seu,Johnson2021.seu,Abdelfattah2022.seu))
x10.seu$platform <- "10x"
x10.seu <- JoinLayers(x10.seu)
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
x10.seu <- standard10X(x10.seu,verbose = T,res = 0.2)
aa <- FindAllMarkers(x10.seu,logfc.threshold = 2,only.pos = T)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DimPlot(x10.seu,group.by = "orig.ident",label = T)+DimPlot(x10.seu,label = T)

qs::qsave(x10.seu,file = "10x_int.seu.sq")
qs::qsave(ss2.seu,file = "ss2_int.seu.sq")

int.seu <- merge(x10.seu,ss2.seu)
int.seu <- JoinLayers(int.seu)
int.seu <- standard10X(int.seu,verbose = T)
