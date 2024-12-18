##Subcluster##
.libPaths()
setwd("~/rna/sc/tarun")
rm(list=ls())
gc()
library(Seurat)
library(dplyr)
packageVersion("Seurat")
path <- dir("./raw/")
sample_names <- path
scRNAlist <- list()
for (i in 1:length(path)) {
  tmp <- Read10X(paste0("./raw/",path[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = tmp,
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)
int.seu <- subset(int.seu, subset = nFeature_RNA >=3000&percent.mt<10&nCount_RNA <150000)
FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

int.seu = JoinLayers(int.seu)
int.seu = NormalizeData(int.seu,verbose=T,assay = "RNA")
int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=T)
int.seu = ScaleData(int.seu,verbose=T)
int.seu = RunPCA(int.seu,verbose=T)
ElbowPlot(int.seu,ndims = 50)
int.seu <- harmony::RunHarmony(int.seu, group.by.vars = "orig.ident",
                               reduction = "pca", assay.use = "RNA",
                               reduction.save = "harmony",max_iter=20)
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:15)
int.seu = FindNeighbors(int.seu,dims=seq(25),verbose=F,reduction = "pca")
int.seu <- FindClusters(int.seu,resolution = 0.3)
VlnPlot(int.seu,features = c("nFeature_RNA","nCount_RNA"))
DimPlot(int.seu,label = T)+DimPlot(int.seu,group.by = "orig.ident",label = T)
FeaturePlot(int.seu,features = c("EGFR","MET","PDGFRA","CDK4","STK17A"),order = T)&scale_colour_viridis_c(option = "E")

all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 1) 

markers <- FindMarkers(int.seu,ident.1 = "131_c1",ident.2 = "131_c2")
marker2 <- FindMarkers(int.seu,ident.1 = "1215",ident.2 = "1215_c2")

top10.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

int.seu <- Seurat::CellCycleScoring(int.seu,s.features = cc.genes.updated.2019$s.genes,
                                    g2m.features = cc.genes.updated.2019$g2m.genes)
DimPlot(int.seu,group.by = "Phase")

qs::qsave(int.seu,"int.seu.qs")
#######################
int.seu <- qs::qread("int.seu.qs")

seurat_colors <- sample(MetBrewer::met.brewer("Klimt",23),23)
features <- list(
  mo  = c("NF2","ADGRG6","SOX10","FOXD3"),
  m1  = c("YAP1","WWTR1","AXL","MERTK","TYRO3","MET"),
  m3  = c("SOX2","NES","CD3G","MS4A1","AIF1","TMEM119"))
DimPlot(int.seu,reduction = "umap",
        label = T,repel = T,label.size = 4,label.box = T,
        cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)+
  DotPlot(int.seu,features = features)
DotPlot(int.seu,features = "VWF")
library(ggplot2)
FeaturePlot(int.seu,features = unlist(features),order = T,
            reduction = "umap",raster = F)&scale_colour_viridis_c(option = "B")
ggsave("vst5_featureplot.png",height = 12,width = 16)

features <-c("EPCAM","PECAM1","MME","NRXN1","ITGB8","CD3G","CD3E","CD79A","PTPRC")
FeaturePlot(int.seu,features = unlist(features),order = T,
            reduction = "umap",raster = F)&scale_colour_viridis_c(option = "B")
ggsave("vst5_featureplot2.png",height = 12,width = 16)
int.seu <- JoinLayers(int.seu)


cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker.csv")

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",10),10)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 3,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP.png",dpi = 300,width = 8,height = 6)
DimPlot(int.seu, group.by = "orig.ident",reduction = "umap", label = T,repel = T,label.size = 2,label.box = T,cols = sample(MetBrewer::met.brewer("Klimt",4),4),label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP_subtype.png",dpi = 300,width = 12,height = 6)
qs::qsave(int.seu,file = "./int.seu.qs")

DimPlot(int.seu,group.by = c("seurat_clusters","celltype"))

int.seu <- qs::qread("int.seu.qs")


DimPlot(int.seu)

# GBM1215 <- subset(int.seu,orig.ident=="1215")
# standard10X <- function(int.seu,pc=30,res=0.3){
# int.seu = NormalizeData(int.seu,verbose=T,assay = "RNA")
# int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=T)
# int.seu = ScaleData(int.seu,verbose=T)
# int.seu = RunPCA(int.seu,verbose=T)
# ElbowPlot(int.seu,ndims = 50)
# # int.seu <- harmony::RunHarmony(int.seu, group.by.vars = "orig.ident", 
# #                                reduction = "pca", assay.use = "RNA", 
# #                                reduction.save = "harmony",max_iter=20)
# int.seu <- RunUMAP(int.seu,reduction = "pca",dims = seq(pc))
# int.seu = FindNeighbors(int.seu,dims=seq(pc),verbose=F,reduction = "pca")
# int.seu <- FindClusters(int.seu,resolution = res)
# }
# GBM1215 <- standard10X(GBM1215,pc=10)
# GBM914 <- standard10X(subset(int.seu,orig.ident=="914"))
# GBM131 <- standard10X(subset(int.seu,orig.ident=="131"))
# GBM1219 <- standard10X(subset(int.seu,orig.ident=="1219"))
# 
# DimPlot(GBM1215)
# all.markers <- FindAllMarkers(GBM1215, only.pos = TRUE, min.pct = 0.25,
#                               logfc.threshold = 1) 
