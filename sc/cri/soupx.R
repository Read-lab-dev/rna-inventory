##Subcluster##
setwd("/home/hzg/rna2/cri")
rm(list=ls())
gc()
int.seu1 <-qs::qread("tumor.qs")
#####
library(SoupX)
library(Seurat)

sc1 = autoEstCont(load10X("~/rna2/091423sc/21047-120-01-01-01/outs"))
sc2 = autoEstCont(load10X("~/rna2/091423sc/21047-120-02-01-01/outs"))
sc3 = autoEstCont(load10X("~/rna2/091423sc/21047-120-03-01-01/outs"))
out1 = adjustCounts(sc1)
out2 = adjustCounts(sc2)
out3 = adjustCounts(sc3)
scRNAlist <- list()
adjcount_names <- c("out1","out2","out3")
sample_names <- c("CrizL","CrizS","DMSO")
for (i in 1:3) {
  scRNAlist[[i]] <- CreateSeuratObject(counts = adjcount_names[i],
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
out1 <- CreateSeuratObject(counts = out1,project = sample_names[1],min.cells=3,min.features = 200)
out2 <- CreateSeuratObject(counts = out2,project = sample_names[2],min.cells=3,min.features = 200)
out3 <- CreateSeuratObject(counts = out3,project = sample_names[3],min.cells=3,min.features = 200)

int.seu <- merge(out1,c(out2,out3))

int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)

int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.85&percent.mt<20& 
                    nFeature_RNA >=250&nCount_RNA>=1000)

qs::qsave(int.seu,file = "soupx.qs")

standard10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  int.seu = NormalizeData(dat,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 2000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose) %>% identity()
  return(int.seu)
}

int.seu <- standard10X(int.seu, nPCs=20, res=0.6)
Idents(int.seu) <- int.seu$seurat_clusters
DimPlot(int.seu)
qs::qsave(int.seu,file = "int_adj.qs")

org.seu <- qs::qread("org.seu.qs")
gsc.seu <- qs::qread("gsc.seu.qs")
