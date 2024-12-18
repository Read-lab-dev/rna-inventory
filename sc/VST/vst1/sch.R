##Subcluster##
.libPaths()
setwd("~/rna/sc/VST/vst1/")
rm(list=ls())
gc()
library(Seurat)
library(dplyr)
packageVersion("Seurat")
path <- dir("./raw/")
sample_names <- gsub(".*\\_","",path)
scRNAlist <- list()
for (i in 1:length(path)) {
  tmp <- Read10X(paste0("./raw/",path[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = tmp,
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
meta <- read.csv("scMeta.csv",row.names = 1)
table(rownames(meta)%in%colnames(int.seu))
int.seu <- subset(int.seu,cells = rownames(meta))
int.seu@meta.data <- meta
int.seu = NormalizeData(int.seu,verbose=F)
int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=F)
int.seu = ScaleData(int.seu,verbose=F)
int.seu = RunPCA(int.seu,verbose=F)
ElbowPlot(int.seu,ndims = 50)
int.seu <- harmony::RunHarmony(int.seu, group.by.vars = c("orig.ident"), 
                               reduction = "pca", assay.use = "RNA", 
                               reduction.save = "harmony",max_iter=20)
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu <- RunTSNE(int.seu,reduction = "harmony",dims = 1:30)
int.seu = FindNeighbors(int.seu,dims=seq(30),verbose=F,reduction = "harmony")
Idents(int.seu) <- int.seu$final_label
#######################
int.seu <- qs::qread("vst1.qs")
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",12),12)
DimPlot(int.seu, group.by = "final_label",reduction = "umap", 
        label = T,repel = T,label.size = 4,label.box = T,
        cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)

DimPlot(int.seu, group.by = "technique",reduction = "umap", 
        label = T,repel = T,label.size = 4,label.box = T,
        cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)


features <- list(
  mo  = c("NF2","ADGRG6","SOX10","FOXD3"),
  m1  = c("YAP1","WWTR1","AXL","MERTK","TYRO3","MET"),
  m3  = c("SOX2","NES","CD3G","MS4A1","AIF1","TMEM119"))

library(ggplot2)
FeaturePlot(int.seu,features = unlist(features),order = T,
            reduction = "umap",raster = F)&scale_colour_viridis_c(option = "B")

DotPlot(int.seu,features = features)

DimPlot(int.seu, group.by = "final_label",reduction = "umap", 
        label = T,repel = T,label.size = 4,label.box = T,
        cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)

ggsave("UMAP.png",dpi = 300,width = 8,height = 6)

DimPlot(int.seu, group.by = "technique",reduction = "umap", 
        label = T,repel = T,label.size = 4,label.box = T,
        cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)

ggsave("UMAP_tech.png",dpi = 300, width = 8,height = 6)

features <- list(
  mo  = c("NF2","ADGRG6","SOX10","FOXD3"),
  m1  = c("YAP1","WWTR1","AXL","MERTK","TYRO3","MET"),
  m3  = c("SOX2","NES","CD3G","MS4A1","AIF1","TMEM119"))
FeaturePlot(int.seu,features = unlist(features),order = T,
            reduction = "umap",raster = F)&scale_colour_viridis_c(option = "B")
ggsave("vst1_featureplot.png",height = 12,width = 16)


features <-c("EPCAM","PECAM1","MME","NRXN1","ITGB8","CD3G","CD3E","CD79A","PTPRC")
FeaturePlot(int.seu,features = unlist(features),order = T,
            reduction = "umap",raster = F)&scale_colour_viridis_c(option = "B")
ggsave("vst1_featureplot2.png",height = 12,width = 16)

DimPlot(int.seu, group.by = "chr22q_loss",reduction = "umap", 
        label = T,repel = T,label.size = 4,label.box = T,
        cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)


qs::qsave(int.seu,file = "./sch.seu.qs")
