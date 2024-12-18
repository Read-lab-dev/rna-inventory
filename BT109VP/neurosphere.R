##02-08-2024##
##Neurosphere##
setwd("~/backup/lab/BT109VP")
rm(list=ls())
gc()
library(SoupX)
library(Seurat)
library(dplyr)
path <- grep(pattern="21047FL-89-01", dir("/home/hzg/backup/lab/ProcessedBam/BT109_VP"),value = T)
sample_names <- c("BT109DMSO","BT109VP")
out1 <- Seurat::Read10X("/home/hzg/backup/21047FL-89/21047FL-89-01-05/outs/filtered_feature_bc_matrix/")
out2 <- Seurat::Read10X("/home/hzg/backup/21047FL-89/21047FL-89-01-06/outs/filtered_feature_bc_matrix/")

scRNAlist <- list()
for (i in 1:2) {
  scRNAlist[[i]] <- CreateSeuratObject(counts = get(paste0("out",i)),
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)

int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)

int.seu <- subset(int.seu, subset = percent.mt<25& 
                    nFeature_RNA >=250&nCount_RNA>=1000)
##Normal###
standard10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  int.seu = NormalizeData(dat,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}

int.seu <-standard10X(int.seu)
DimPlot(int.seu,group.by = "orig.ident")
FeaturePlot(int.seu,features = "H3K27M")
###SCT##
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
int.seu <- CellCycleScoring(int.seu,s.features = s.genes,g2m.features = g2m.genes)
int.seu <- SCTransform(int.seu,vars.to.regress =  c("percent.mt","S.Score","G2M.Score"),verbose = F)
DefaultAssay(int.seu) <-"SCT"
int.seu <- RunPCA(int.seu,assay="SCT")
ElbowPlot(int.seu)
int.seu <- RunUMAP(int.seu,reduction = "pca",dims = 1:30)
int.seu <- RunTSNE(int.seu,reduction = "pca",dims = 1:30)
int.seu <- FindNeighbors(int.seu,dims=seq(30),verbose=F,reduction = "pca")
int.seu <- FindClusters(int.seu,res=0.5,verbose=F)
DimPlot(int.seu,reduction = "tsne")
all.markers <-FindMarkers(int.seu,ident.1 = "BT109VP",ident.2 = "BT109DMSO",group.by="orig.ident")
all.markers <-FindAllMarkers(int.seu,logfc.threshold = 0.5)

top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
qs::qsave(int.seu,"neurosphere.qs")

neu.seu <- qs::qread("neurosphere.qs")


path <- grep(pattern="21047FL-89-01", dir("/home/hzg/backup/lab/ProcessedBam/BT109_VP"),value = T)
sample_names <- c("BT109DMSO","BT109VP")
out1 <- Seurat::Read10X("/home/hzg/backup/21047FL-89/21047FL-89-01-05/outs/filtered_feature_bc_matrix/")
out2 <- Seurat::Read10X("/home/hzg/backup/21047FL-89/21047FL-89-01-06/outs/filtered_feature_bc_matrix/")

scRNAlist <- list()
for (i in 1:2) {
  scRNAlist[[i]] <- CreateSeuratObject(counts = get(paste0("out",i)),
                                       project = sample_names[i])
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
int.seu <- int.seu %>% NormalizeData() %>% ScaleData()

mut <- as.data.frame(t(as.data.frame(int.seu@assays[["RNA"]]@counts)))

mut_merge <- left_join(tibble::rownames_to_column(neu.seu@meta.data), tibble::rownames_to_column(mut))

neu.seu$H3K27M <- mut_merge$H3K27M
neu.seu$H3K27wt <- mut_merge$H3K27WT
neu.seu$H3G34R <- mut_merge$H3G34R

FeaturePlot(neu.seu,features = c("H3K27M","H3K27wt","H3G34R"),split.by = "orig.ident")


#########Detect H3K27M########
int.seu$k27 <- ifelse(GetAssayData(int.seu)["H3K27M",]>0,'K27M','K27')
int.seu$barcode <- colnames(int.seu)
int.seu$gfp <- ifelse(GetAssayData(int.seu)["eGFP",]>0,'eGFP+','eGFP-')
table(int.seu$k27,int.seu$gfp)
DimPlot(int.seu,cells.highlight = intersect(colnames(int.seu),list$V1))
FeaturePlot(int.seu,features = c("H3K27M","eGFP"),split.by = "orig.ident")
