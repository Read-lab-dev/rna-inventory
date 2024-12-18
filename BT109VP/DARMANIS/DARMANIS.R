##
rm(list = ls())
setwd("~/rna/BT109VP/DARMANIS")
library(tibble)
exp <- data.table::fread("GBM_raw_gene_counts.csv.gz")
exp <- as.data.frame(exp)
exp <- column_to_rownames(exp,var = "V1")
gbm.seu <- CreateSeuratObject(counts = exp)
meta.data <- read.table("GBM_metadata.csv")
gbm.seu@meta.data <- cbind(gbm.seu@meta.data,meta.data)

###SCT##
gbm.seu[["percent.mt"]] <- PercentageFeatureSet(gbm.seu, pattern = "^MT-")
gbm.seu$orig.ident <- gbm.seu$sample
library(harmony)
standard10X = function(int.seu,nPCs=30,res=0.8,verbose=FALSE){
  int.seu = NormalizeData(int.seu,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures=3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu,group.by.vars="sample", max_iter = 20)
  int.seu = RunUMAP(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
gbm.seu <- standard10X(gbm.seu)
# tsne <- read.table("GBM_TSNE.csv")
# colnames(tsne) <- c("tSNE_1","tSNE_2")
# gbm.seu@reductions$tsne@cell.embeddings <- as.matrix(tsne)
# DimPlot(gbm.seu,group.by = "Selection",reduction = "tsne")
# DimPlot(gbm.seu,reduction = "tsne")
marker <- FindAllMarkers(gbm.seu,logfc.threshold = 1)
top10.markers <- marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
Idents(gbm.seu) <- gbm.seu$celltype
DimPlot(gbm.seu,label = T,label.box = T)&NoLegend()
saveRDS(gbm.seu,file = "darmanis.seu.Rds")


