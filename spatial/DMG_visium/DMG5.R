rm(list = ls())
.libPaths(c("/home/hzg/R/x86_64-pc-linux-gnu-library/4.2/SeuratV5",.libPaths()))
options(future.globals.maxSize = 8000 * 1024^2)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
packageVersion("Seurat")
setwd("/home/hzg/rna/visium/Org/")
vis.seu <- Load10X_Spatial(data.dir = "./data/DMG5_spaceranger_out",slice ="DMG5") 
vis.seu$orig.ident <- "DMG5"
vis.seu <- SCTransform(vis.seu, assay = "Spatial", verbose = FALSE)

DefaultAssay(vis.seu) <- "SCT"
VariableFeatures(vis.seu)
vis.seu <- RunPCA(vis.seu,assay = "SCT", verbose = FALSE)
vis.seu <- RunUMAP(vis.seu, reduction = "pca",dims = 1:10)
vis.seu <- FindNeighbors(vis.seu,reduction = "pca", dims = 1:10)
vis.seu <- FindClusters(vis.seu, verbose = FALSE)

DimPlot(vis.seu, reduction = "umap", group.by = c("ident", "orig.ident"))
SpatialDimPlot(vis.seu,group.by = "seurat_clusters")&NoLegend()

# Idents(vis.seu) <- vis.seu$orig.ident
# plot1 <- VlnPlot(vis.seu, features = "nCount_Spatial", pt.size = 0.1) &NoLegend()
# plot2 <- SpatialFeaturePlot(vis.seu, features = "nCount_Spatial")&theme(legend.position = "right")
# wrap_plots(plot1, plot2)

vis.seu <- SCTransform(vis.seu, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(vis.seu, features = c("SOX4", "SOX2"))

vis.seu <- PrepSCTFindMarkers(vis.seu)
de_markers <- FindAllMarkers(vis.seu,logfc.threshold=1)

vis.seu <- FindSpatiallyVariableFeatures(vis.seu,assay = "SCT", features = VariableFeatures(vis.seu)[1:1000],
                                         selection.method = "moransi")

top.features = rownames(
  dplyr::slice_min(
    vis.seu[["SCT"]]@meta.features,
    moransi.spatially.variable.rank,n = 6))
SpatialFeaturePlot(object = vis.seu, features = top.features, ncol = 3)

qs::qsave(vis.seu,file = "DMG5.qs")

###Intergrate:

allen_reference <- qs::qread("~/rna/BT109VP/gsc.seu.qs")
# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
allen_reference <- UpdateSeuratObject(allen_reference)
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
# After subsetting, we renormalize cortex
DimPlot(allen_reference, group.by = "celltype", label = TRUE,reduction = "dim2")

anchors <- FindTransferAnchors(reference = allen_reference, query = vis.seu, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$celltype, prediction.assay = TRUE,
                                  weight.reduction = vis.seu[["pca"]], dims = 1:30)
vis.seu[["predictions"]] <- predictions.assay
DefaultAssay(vis.seu) <- "predictions"

SpatialFeaturePlot(vis.seu, features =  c("tOPC", "tCYC","tNPC","tAC"), pt.size.factor = 1.8, ncol = 2, crop = T, alpha = c(0.1, 1))
