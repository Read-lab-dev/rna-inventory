##########Spitzer2024#########
rm(list = ls())
library(Seurat)
library(dplyr)
setwd("~/rna/sc/LGG")
exp <- data.table::fread("./Spitzer2024/GSE260928_Chromium_Astrocytoma1_processed_UMI.tsv.gz")
exp <- as.data.frame(exp)
rownames(exp)<- exp$V1
exp <- exp[,-1]

exp2 <- data.table::fread("./Spitzer2024/GSE260928_Chromium_Astrocytoma2_processed_UMI.tsv.gz")
exp2 <- as.data.frame(exp2)
rownames(exp2)<- exp2$V1
exp2 <- exp2[,-1]

exp3 <- data.table::fread("./Spitzer2024/GSE260928_Chromium_Oligodendroglioma_processed_UMI.tsv.gz")
exp3 <- as.data.frame(exp3)
rownames(exp3)<- exp3$V1
exp3 <- exp3[,-1]

identical(rownames(exp2),rownames(exp))

A1.seu <- CreateSeuratObject(counts = exp)
A2.seu <- CreateSeuratObject(counts = exp2)
O1.seu <- CreateSeuratObject(counts = exp3)

Spitzer2024.seu <- merge(A1.seu,y=c(A2.seu,O1.seu))
Spitzer2024.seu <- JoinLayers(Spitzer2024.seu)
Spitzer2024.seu$treat <- stringr::str_split(colnames(Spitzer2024.seu),'-',simplify = T)[,1]
Spitzer2024.seu$treat <- stringr::str_split(Spitzer2024.seu$treat,'_',simplify = T)[,2]
Spitzer2024.seu$patient <- Spitzer2024.seu$orig.ident
Spitzer2024.seu$orig.ident <- stringr::str_split(colnames(Spitzer2024.seu),'-',simplify = T)[,1]
table(Spitzer2024.seu$orig.ident)

Spitzer2024.seu[["percent.mt"]] <- PercentageFeatureSet(Spitzer2024.seu, pattern = "^MT-")
VlnPlot(Spitzer2024.seu, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
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
Spitzer2024.seu <- standard10X(Spitzer2024.seu,verbose = T)
DimPlot(Spitzer2024.seu,group.by = "orig.ident",label = T)+DimPlot(Spitzer2024.seu,label = T)
aa <- FindAllMarkers(Spitzer2024.seu,logfc.threshold = 2,only.pos = T)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

qs::qsave(Spitzer2024.seu,file = "Spitzer2024-10x.seu.qs")
