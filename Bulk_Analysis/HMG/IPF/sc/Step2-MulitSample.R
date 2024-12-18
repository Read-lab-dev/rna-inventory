rm(list = ls())

setwd("/home/hzg/rna/scRNA_analysis/IPF/sc")
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)

##Load the dataset you want to merge

path <- grep(pattern="h5", dir(),value = T)

sample_names <- gsub("_filtered_gene_bc_matrices_h5.h5","",path)

sample_names <- gsub("GSM......._","",sample_names)
# Load and QC
scRNAlist <- list()
for (i in 1:length(path)) {
  testdf <- Read10X_h5(path[i])
  
  scRNAlist[[i]] <- CreateSeuratObject(counts = testdf,project = sample_names[i],min.cells=3,min.features = 200)
  
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id= sample_names[i])
  
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
}
# Simple Merge
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)

FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")+
  FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.85&percent.mt<20& 
                    nFeature_RNA >=250&nCount_RNA>=1000)
# Standard PCA
int.seu <- int.seu %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst",nfLeatures= 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50,verbose = F)

ElbowPlot(int.seu, ndims=50, reduction="pca") 
pc.num <- 1:20

int.seu <- int.seu %>% 
  RunUMAP(reduction = "pca", dims = pc.num, verbose = F) %>% 
  FindNeighbors(reduction = "pca", k.param = 10, dims = pc.num)
DimPlot(int.seu,group.by = "type")

int.seu$type <- gsub("_..","",int.seu$orig.ident)

hpca.se <-celldex::HumanPrimaryCellAtlasData()
assay_for_SingleR <- GetAssayData(int.seu, slot="data") 
int.hesc <- SingleR(test = assay_for_SingleR, ref = hpca.se, labels = hpca.se$label.main) 
int.hesc$labels

int.seu$singleR <- int.hesc$labels
Idents(int.seu) <- int.seu$singleR
int.seu <-subset(int.seu,idents=names(table(int.seu$singleR))[table(int.seu$singleR)>50])

DimPlot(int.seu,cols = sample(met.brewer("Klimt",14),14))
  ##Save the merged dataset

saveRDS(int.seu,file="Merged_All_none_integrated.Rds")


