rm(list = ls())
.libPaths(c("/home/hzg/R/x86_64-pc-linux-gnu-library/4.2/SeuratV5",.libPaths()))
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
packageVersion("Seurat")
setwd("/home/hzg/backup/lab/visium/Org/")
list.files("./tar") # Should show filtered_feature_bc_matrix.h5
for (i in 2:6) {
  idx <- list.files("./tar")[i]
  vis.seu <- Load10X_Spatial(data.dir = paste0(getwd(),"/tar/",idx),slice =substr(idx,1,4) ) 
  vis.seu$orig.ident <- substr(idx,1,4)
  vis.seu <- SCTransform(vis.seu, assay = "Spatial", verbose = FALSE)
  assign(substr(idx,1,4),vis.seu)
}
brain.merge <- merge(DMG1, c(DMG2,DMG3,DMG4,DMG5))

DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(DMG1), VariableFeatures(DMG2),VariableFeatures(DMG3),
                                   VariableFeatures(DMG4),VariableFeatures(DMG5))
brain.merge <- RunPCA(brain.merge,assay = "SCT", verbose = FALSE)
library(harmony)
brain.merge <- RunHarmony(brain.merge,group.by.vars="orig.ident")
brain.merge <- RunUMAP(brain.merge, reduction = "harmony",dims = 1:30)
brain.merge <- FindNeighbors(brain.merge,reduction = "harmony", dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)

DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
SpatialDimPlot(brain.merge,group.by = "seurat_clusters",ncol = 3)&NoLegend()

rm(DMG1,DMG2,DMG3,DMG4,DMG5,vis.seu)
gc()

Idents(brain.merge) <- brain.merge$orig.ident
plot1 <- VlnPlot(brain.merge, features = "nCount_Spatial", pt.size = 0.1) &NoLegend()
plot2 <- SpatialFeaturePlot(brain.merge, features = "nCount_Spatial")&theme(legend.position = "right")
wrap_plots(plot1, plot2)
brain.merge <- SCTransform(brain.merge, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(brain.merge, features = c("SOX4", "SOX2"))

p1 <- SpatialFeaturePlot(brain.merge, features = "TTR", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain.merge, features = "TTR", alpha = c(0.1, 1))
p1 + p2

brain.merge <- PrepSCTFindMarkers(brain.merge)
de_markers <- FindAllMarkers(brain.merge,logfc.threshold=1)
top10.markers <- de_markers[-c(grep(c("^RP[SL]"),rownames(de_markers)),
                                grep(c("^MT-"),rownames(de_markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
SpatialFeaturePlot(object = brain.merge, features = c("SOX2","MBP"), alpha = c(0.1, 1), ncol = 2)

brain.merge <- FindSpatiallyVariableFeatures(brain.merge,assay = "SCT", features = VariableFeatures(brain.merge)[1:1000],
                                       selection.method = "moransi")

top.features <- head(SpatiallyVariableFeatures(brain.merge, selection.method = "moransi"), 20)
SpatialFeaturePlot(brain.merge, features = c("OLIG1","GFAP","RBFOX3"), alpha = c(0.1, 1))
qs::qsave(brain.merge,file = "brain.merge.qs")

cortex <- subset(brain.merge, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
| image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2

###Intergrate:

allen_reference <- readRDS("/brahms/shared/vignette-data/allen_cortex.rds")
# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "moransi",
                                        features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex, selection.method = "moransi"), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)

SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT",
                                        "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))