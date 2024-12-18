rm(list=ls())
gc()
library(Seurat)
library(dplyr)
packageVersion("Seurat")
path <- dir("./raw/")
sample_names <- path
scRNAlist <- list()
for (i in 1:length(path)) {
  tmp <- Read10X(paste0("./raw/",path[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = tmp,
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
int.seu <- JoinLayers(int.seu)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)
int.seu <- subset(int.seu, subset = nFeature_RNA >=250&percent.mt<20)
int.seu = NormalizeData(int.seu,verbose=F,assay = "RNA")
int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=F)
int.seu = ScaleData(int.seu,verbose=F)
int.seu = RunPCA(int.seu,verbose=F)
ElbowPlot(int.seu,ndims = 50)
int.seu <- harmony::RunHarmony(int.seu, group.by.vars = "orig.ident", 
                               reduction = "pca", assay.use = "RNA", 
                               reduction.save = "harmony",max_iter=20)
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu = FindNeighbors(int.seu,dims=seq(30),verbose=F,reduction = "harmony")
int.seu <- FindClusters(int.seu,resolution = 0.4)
DimPlot(int.seu)
saveRDS(int.seu,file = "Adam.Rdata")

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SingleCellExperiment)
reducedDim(SCPCP000001_merged) <- NULL
counts(sce)
manno.seurat <- as.Seurat(SCPCP000001_merged, counts = "counts")

