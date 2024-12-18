library(NMF)
gsc.seu <- qs::qread("int.seu.qs")
gsc.seu <- subset(gsc.seu,seurat_clusters==9,invert=T)
gsc.seu <- subset(gsc.seu, downsample = 800)
gsc.seu = NormalizeData(gsc.seu,verbose=F)
gsc.seu = FindVariableFeatures(gsc.seu,selection.method = "vst",nfeatures= 3000,verbose=F)
gsc.seu = ScaleData(gsc.seu,verbose=F)
mat <- gsc.seu@assays$RNA@scale.data
mat <- as.matrix(gsc.seu@assays$RNA@counts[rownames(gsc.seu@assays$RNA@scale.data),])
mat <- mat[rowSums(mat)>0,]
res.rank <- NMF::nmf(mat,
                     rank = 4:10,
                     .opt="vP10",
                     method = "snmf/r")
plot(res.rank)

nmf.rank6 <- NMF::nmf(mat,
                      rank = 6,
                      nrun=50,
                      .opt="vP10",
                      method = "snmf/r")
save(gsc.seu,res.rank,file = "res.rank.Rdata")
.libPaths(c("/home/hzg/R/x86_64-pc-linux-gnu-library/SeuratV5/",.libPaths()))
library(singlet)
library(dplyr)
library(cowplot)
set.seed(123) # for reproducible NMF models
gsc.seu <- NormalizeData(gsc.seu)
gsc.seu <- RunNMF(gsc.seu,k=6)
gsc.seu <- RunUMAP(gsc.seu, reduction = "nmf", dims = 1:ncol(gsc.seu@reductions$nmf))
gsc.seu <- FindNeighbors(gsc.seu, dims = 1:ncol(gsc.seu@reductions$nmf), reduction = "nmf") %>%
  FindClusters(resolution = 0.5, verbose = FALSE)
DimPlot(gsc.seu,reduction = "umap",label = T)
gsc.seu <- calculate_state(gsc.seu)
DimPlot(gsc.seu,reduction = "dim2")
DimPlot(gsc.seu,group.by = "state")
