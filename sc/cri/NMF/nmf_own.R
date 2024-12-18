devtools::install_github("zdebruine/RcppML")
devtools::install_github("zdebruine/singlet")
.libPaths(c("/home/hzg/R/x86_64-pc-linux-gnu-library/SeuratV5/"))
library(RcppML)
library(Seurat)
library(singlet)
int.seu <- qs::qread("int.seu.qs")
int.seu <- subset(int.seu,seurat_clusters==9,invert=T)
library(dplyr)
library(cowplot)
set.seed(123) # for reproducible NMF models
int.seu <- NormalizeData(int.seu)
int.seu <- FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000)
int.seu <- int.seu[rownames(int.seu@assays$RNA$scale.data),]
int.seu <- RunNMF(int.seu,k=8)
int.seu <- RunUMAP(int.seu, reduction = "nmf", dims = 1:ncol(int.seu@reductions$nmf))
int.seu <- FindNeighbors(int.seu, dims = 1:ncol(int.seu@reductions$nmf), reduction = "nmf") %>%
  FindClusters(resolution = 0.5, verbose = FALSE)
singlet::RankPlot(int.seu)
DimPlot(int.seu,reduction = "umap",label = T,group.by = "state")

DimPlot(int.seu,reduction = "umap",label = T,group.by = "seurat_clusters")
DimPlot(int.seu,reduction = "umap",label = T,group.by = "state")+
  DimPlot(int.seu,reduction = "umap",label = T,group.by = "nmf_cluster")
MetadataPlot(int.seu,"state",reduction = "nmf")+DimPlot(int.seu,group.by = c("state"),label = T)
lineage <- read.csv("./own6.csv",header = T,na.strings = "")[1:50,]
calculate_state <- function(int.seu){
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
  library(ggsci)
  library(gghighlight)
  library(MetBrewer)
  library(tibble)
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$MES1)),nbin = 30,ctrl = 100,name = "MES1")
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$MES2)),nbin = 30,ctrl = 100,name = "MES2")
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$AC)),nbin = 30,ctrl = 100,name = "AC")
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$OPC)),nbin = 30,ctrl = 100,name = "OPC")
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$NPC)),nbin = 30,ctrl = 100,name = "NPC1")
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$NPC)),nbin = 30,ctrl = 100,name = "NPC2")
  
  gsc.module <- data.frame(int.seu$MES11,int.seu$MES21,int.seu$AC1,
                           int.seu$OPC1,int.seu$NPC11,int.seu$NPC21)
  colnames(gsc.module) <- c("MES1","MES2","AC","OPC","NPC1","NPC2")
  
  gsc.module$state <- apply(gsc.module,1,which.max)
  gsc.module$state <- gsub(1,"MES",gsc.module$state)
  gsc.module$state <- gsub(2,"MES",gsc.module$state)
  gsc.module$state <- gsub(3,"AC",gsc.module$state)
  gsc.module$state <- gsub(4,"OPC",gsc.module$state)
  gsc.module$state <- gsub(5,"NPC",gsc.module$state)
  gsc.module$state <- gsub(6,"NPC",gsc.module$state)
  gsc.module$MES <- apply(gsc.module[,1:2],1,max)
  gsc.module$NPC <- apply(gsc.module[,5:6],1,max)
  gsc.module <-gsc.module[,c(8,9,3,4,7)]
  gsc.module$D <- ifelse(gsc.module$state%in%c("NPC","MES"),0,1)
  yaxis <- function(x){(max(x[4],x[2])-max(x[3],x[1]))}
  
  xaxis <- function(x){aa =ifelse(x[5]>0, log2(abs(x[4]-x[2])+1),
                                  log2(abs(x[3]-x[1])+1))
  aa= ifelse(x[6]>0,-aa,aa)
  return(aa)
  }
  gsc.module$yaxis <- apply(gsc.module[,c(1:4)],1,yaxis)
  gsc.module$xaxis <- apply(gsc.module[,c(1:4,7,6)], 1, xaxis)
  gsc.module$ident <- int.seu$orig.ident
  int.seu$state <- gsc.module$state
  redc <- as.matrix(gsc.module[,8:7])
  colnames(redc) <- c("Dim2_1","Dim2_2")
  int.seu[["dim2"]] <- CreateDimReducObject(embeddings = redc,key = "Dim2_")
  DimPlot(int.seu,reduction = "dim2")
  return(int.seu)
}
int.seu <- calculate_state(int.seu)
DimPlot(int.seu,reduction = "dim2",group.by = "state")
