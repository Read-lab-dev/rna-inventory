rm(list = ls())

.libPaths(c("/home/hzg/R/x86_64-pc-linux-gnu-library/SeuratV5/",.libPaths()))
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
packageVersion("Seurat")
setwd("~/rna/visium")
vis.seu <- Load10X_Spatial(data.dir =getwd())
vis.seu$orig.ident <- "vis.seu"
vis.seu <- SCTransform(vis.seu, assay = "Spatial", verbose = FALSE)
DefaultAssay(vis.seu) <- "SCT"
# vis.seu <- qs::qread("vis.seu.qs")
plot1 <- VlnPlot(vis.seu, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(vis.seu, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

SpatialFeaturePlot(vis.seu, features = c("SOX2","OLIG1"))

vis.seu <- RunPCA(vis.seu, assay = "SCT", verbose = FALSE)
vis.seu <- FindNeighbors(vis.seu, reduction = "pca", dims = 1:30)
vis.seu <- FindClusters(vis.seu, verbose = FALSE)
vis.seu <- RunUMAP(vis.seu, reduction = "pca", dims = 1:30)

p1 <- DimPlot(vis.seu, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(vis.seu, label = TRUE, label.size = 3,alpha = 0.5)
p1 + p2

de_markers <- FindAllMarkers(vis.seu,logfc.threshold=1)
top10.markers <- de_markers[-c(grep(c("^RP[SL]"),rownames(de_markers)),
                               grep(c("^MT-"),rownames(de_markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
SpatialFeaturePlot(object = vis.seu, features = c("SOX2","MBP"), alpha = c(0.1, 1), ncol = 2)

vis.seu <- FindSpatiallyVariableFeatures(vis.seu, assay = "SCT", features = VariableFeatures(vis.seu)[1:1000],
                                       selection.method = "moransi")

top.features <- head(SpatiallyVariableFeatures(vis.seu, selection.method = "moransi"), 9)
SpatialFeaturePlot(vis.seu, features = top.features, alpha = c(0.1, 1))

setwd("/home/hzg/backup/lab/visium/Org/")
qs::qsave(vis.seu,file = "vis.seu.qs")

##########################################
allen_reference <- qs::qread("/home/hzg/rna/BT109VP/gsc.seu.qs")
allen_reference <- UpdateSeuratObject(allen_reference)
# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
allen_reference <- SCTransform(allen_reference, vars.to.regress = c("orig.ident","percent.mt") ,ncells = 3000, verbose = FALSE)
#   
allen_reference <- NormalizeData(allen_reference) %>% FindVariableFeatures(method = "vst",nfeatures = 3000)
allen_reference<- RunPCA(allen_reference,verbose = FALSE) %>%
  RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>%
  FindClusters(method="igraph",resolution = 0.8,algorithm = 4)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, label = TRUE,reduction = "tsne")
vis.seu <- qs::qread("./seurat/DMG3.qs")

anchors <- FindTransferAnchors(reference = allen_reference, query = vis.seu, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = Idents(allen_reference), prediction.assay = TRUE,
                                  weight.reduction = vis.seu[["pca"]], dims = 1:30)
vis.seu[["predictions"]] <- predictions.assay
DefaultAssay(vis.seu) <- "predictions"
SpatialFeaturePlot(vis.seu, features = rownames(vis.seu),pt.size.factor = 1.6, ncol = 5, crop = TRUE,alpha = c(0.4, 1))
ggsave("projection_Neu_DMG3.png",height = 9,width = 16)
DefaultAssay(vis.seu) <- "SCT"
SpatialFeaturePlot(vis.seu,features = c("OLIG1","EGFR","PROM1","NES","PDGFRA","VIM","SOX2","PAX6","HIF1A","CXCL8"),ncol = 5,pt.size.factor = 1.6,alpha =c(0.1, 1) )
ggsave("FEARTURE2_DMG3.PNG",height = 9,width = 16)

######
# 0         1         2         3         4         5         6         7         8         9        10        11 
# "#544799" "#D2C564" "#406E89" "#D77CA6" "#924099" "#CF5A79" "#DF9ED4" "#D77B5A" "#98B46C" "#3C5096" "#E2AD5F" "#438B7D" 
# 12        13        14 
# "#5DA373" "#CB4955" "#734399" 
SpatialDimPlot(vis.seu,cols = seurat_colors,label = T,label.size = 2,ncol = 3)+
  DimPlot(vis.seu,cols = seurat_colors,label = T,label.box = T,label.size = 2)&NoLegend()

ggsave("Intergrate_DMG_DimPlot.png")
SpatialFeaturePlot(vis.seu,ncol = 3,features = "MIF")+
  FeaturePlot(vis.seu,features = "MIF",order = T)
ggsave("Intergrate_DMG_MIF.png")

SpatialFeaturePlot(vis.seu, features = "4", pt.size.factor = 1.6, ncol = 5, crop = TRUE,alpha = c(0.4, 1),keep.scale = NULL)


########
Tumor
Invasive
Invasive
Invasive
Hypoxic
Tumor
Invasive
Hypoxic
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(vis.seu)
vis.seu <- RenameIdents(vis.seu, new.cluster.ids)
SpatialDimPlot(vis.seu,pt.size.factor = 1.2,label = T)+SpatialFeaturePlot(vis.seu,features = "CHI3L1",pt.size.factor = 1.2)+ggsave("real.png",dpi = 80)

vis.seu <- qs::qread("./Org/seurat/DMG1.qs")
SpatialFeaturePlot(vis.seu, features = rownames(vis.seu)[-5], keep.scale = "feature",pt.size.factor = 1.6, ncol = 4, crop = TRUE,alpha = c(0.4, 1))
ggsave("Intergrate_DMG1_NPC.png")

SpatialFeaturePlot(vis.seu,features="VIM",pt.size.factor = 1.5,alpha=0.8)