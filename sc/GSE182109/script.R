##########Abdelfattah2022#########
rm(list = ls())
library(Seurat)
library(dplyr)
library(reticulate)
library(anndata)
library(Seurat) 
library(Matrix)
setwd("~/rna/sc/GSE182109/")
path <- grep("GSM",dir("./raw"),value = T)
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

FeaturePlot(gbm.seu,reduction = "dim2",split.by = "EGFR_status",features = c("STK17A","CDK4"),order = T,blend = TRUE)
ggsave("STK17A-CDK4-2dim.png",height = 7.5,width = 15)

FeaturePlot(gbm.seu,reduction = "dim2",split.by = "EGFR_status",features = "STK17A")

int.seu <- qs::qread("Abdelfattah2022-10x.seu.qs")
meta <- read.delim2("Meta_Data_GBMatlas.txt")
umap <- read.delim2("UMAP_Data_GBMatlas.txt")
umap$UMAP_1 <- as.numeric(umap$UMAP_1);umap$UMAP_2 <- as.numeric(umap$UMAP_2)
int.seu <- Abdelfattah2022.seu[,rownames(umap)]
umap <- umap[colnames(int.seu),]
meta <- meta[colnames(int.seu),]
library(dplyr)
int.seu[["UMAP"]] <- CreateDimReducObject(embeddings= umap,key = "UMAP_",assay = DefaultAssay(int.seu))
DimPlot(int.seu)
int.seu@meta.data <- meta
DimPlot(int.seu,group.by = "Assignment")
qs::qsave(int.seu,"int.seu.qs")
library(harmony)
mg.seu <- subset(int.seu,Assignment=="Myeloid")
standard10X = function(int.seu,nPCs=30,res=0.3,verbose=FALSE){
  int.seu = NormalizeData(int.seu)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 2000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu,group.by.vars="orig.ident")
  int.seu = RunUMAP(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
mg.seu <- standard10X(mg.seu,res = 0.5)
mg.seu[["UMAP"]] <- NULL
mg.seu <- mg.seu[,!is.na(mg.seu$SubAssignment)]
mg.seu = RunHarmony(mg.seu,group.by.vars=c("Fragment","sex"))
mg.seu <- RunUMAP(mg.seu,reduction = "harmony",dims = 1:10)
DimPlot(mg.seu,group.by = "SubAssignment")
Idents(mg.seu) <- mg.seu$SubAssignment

mg.seu.sub <- subset(mg.seu,downsample=1000)
dat.sub <- as.matrix(mg.seu.sub@assays$RNA$counts)
dat.sub <- dat.sub[rowSums(dat.sub==0)<0.8*ncol(dat.sub),]
info.sub <- t(mg.seu.sub$SubAssignment)
dat.sub <- rbind(info.sub,dat.sub)
write.table(dat.sub,col.names = F,sep = "\t",quote = F)
