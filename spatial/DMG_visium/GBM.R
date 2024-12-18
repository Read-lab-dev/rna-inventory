#GBM Visium:

rm(list = ls())
.libPaths(c("/home/hzg/R/x86_64-pc-linux-gnu-library/4.2/SeuratV5",.libPaths()))
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
packageVersion("Seurat")
setwd("/home/hzg/backup/lab/visium/Org/data/")
list.files(".") # Should show filtered_feature_bc_matrix.h5
for (i in 6:11) {
  idx <- list.files(".")[i]
  vis.seu <- Load10X_Spatial(data.dir = paste0("./",idx),slice =substr(idx,1,4) ) 
  vis.seu$orig.ident <- substr(idx,1,4)
  vis.seu <- SCTransform(vis.seu, assay = "Spatial", verbose = FALSE)
  assign(substr(idx,1,4),vis.seu)
}
brain.merge <- merge(GBM1,c(GBM2,GBM3,GBM4,GBM5,GBM6))

DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(GBM1), VariableFeatures(GBM2),VariableFeatures(GBM3),
                                   VariableFeatures(GBM4),VariableFeatures(GBM5),VariableFeatures(GBM6))
brain.merge <- RunPCA(brain.merge,assay = "SCT", verbose = FALSE)
library(harmony)
brain.merge <- RunHarmony(brain.merge,group.by.vars="orig.ident")
brain.merge <- RunUMAP(brain.merge, reduction = "harmony",dims = 1:30)
brain.merge <- FindNeighbors(brain.merge,reduction = "harmony", dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE,resolution = 0.5)

DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
SpatialDimPlot(brain.merge,group.by = "seurat_clusters",ncol = 3)&NoLegend()

rm(GBM1,GBM2,GBM3,GBM4,GBM5,GBM6,vis.seu)
gc()

Idents(brain.merge) <- brain.merge$orig.ident
plot1 <- VlnPlot(brain.merge, features = "nCount_Spatial", pt.size = 0.1) &NoLegend()
plot2 <- SpatialFeaturePlot(brain.merge, features = "nCount_Spatial")&theme(legend.position = "right")
wrap_plots(plot1, plot2)
brain.merge <- SCTransform(brain.merge, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(brain.merge, features = c("SOX4", "SOX2"))

p1 <- SpatialFeaturePlot(brain.merge, features = "EGFR", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain.merge, features = "EGFR", alpha = c(0.1, 1))
p1 + p2

SpatialFeaturePlot(brain.merge, features = "EGFR",alpha = c(0.1, 1))

slot(object = brain.merge@assays$SCT@SCTModel.list[[6]], name="umi.assay")<-"Spatial"
brain.merge <- PrepSCTFindMarkers(brain.merge)
Idents(brain.merge)<- brain.merge$seurat_clusters
de_markers <- FindAllMarkers(brain.merge,logfc.threshold=1)
top10.markers <- de_markers[-c(grep(c("^RP[SL]"),rownames(de_markers)),
                               grep(c("^MT-"),rownames(de_markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
SpatialFeaturePlot(object = brain.merge, features = c("SOX2","MBP"), alpha = c(0.1, 1), ncol = 6)

brain.merge <- FindSpatiallyVariableFeatures(brain.merge,assay = "SCT", features = VariableFeatures(brain.merge)[1:1000],
                                             selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(brain.merge, selection.method = "moransi"), 20)
SpatialFeaturePlot(brain.merge, features = , alpha = c(0.1, 1))
qs::qsave(brain.merge,file = "GBM.merge.qs")

#########
GBM6 <- RunPCA(GBM6,assay = "SCT", verbose = FALSE)
GBM6 <- RunUMAP(GBM6, reduction = "pca",dims = 1:30)
GBM6 <- FindNeighbors(GBM6,reduction = "pca", dims = 1:30)
GBM6 <- FindClusters(GBM6, verbose = FALSE,resolution = 0.5)

GBM6 <- PrepSCTFindMarkers(GBM6)
GBM6 <- FindSpatiallyVariableFeatures(GBM6,assay = "SCT", features = VariableFeatures(GBM6)[1:1000],
                                             selection.method = "moransi")

top.features <- head(SpatiallyVariableFeatures(GBM6, selection.method = "moransi"), 12)
SpatialFeaturePlot(GBM6, features = top.features, alpha = c(0.1, 1),ncol=6)
ggsave("FEARTURE_GBM6.PNG",height = 9,width = 16)
qs::qsave(GBM6,file = "GBM6.qs")
