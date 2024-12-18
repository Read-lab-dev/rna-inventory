##Subcluster##
.libPaths()
setwd("~/rna/sc/vst3")
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
int.seu <- subset(int.seu, subset = nFeature_RNA >=250&percent.mt<20)
int.seu = NormalizeData(int.seu,verbose=F,assay = "RNA")
int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=F)
int.seu = ScaleData(int.seu,verbose=F)
int.seu = RunPCA(int.seu,verbose=F)
ElbowPlot(int.seu,ndims = 50)
int.seu <- harmony::RunHarmony(int.seu, group.by.vars = "orig.ident", 
                               reduction = "pca", assay.use = "RNA", 
                               reduction.save = "harmony",max_iter=20)
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu = FindNeighbors(int.seu,dims=seq(30),verbose=F,reduction = "harmony")
int.seu <- FindClusters(int.seu,resolution = 0.1)
qs::qsave(int.seu,"sch.seu.qs")

#######################
int.seu <- qs::qread("vst3.qs")
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
ggsave("vst3_featureplot.png",height = 12,width = 16)

features <-c("EPCAM","PECAM1","MME","NRXN1","ITGB8","CD3G","CD3E","CD79A","PTPRC")
FeaturePlot(int.seu,features = unlist(features),order = T,
            reduction = "umap",raster = F)&scale_colour_viridis_c(option = "B")
ggsave("vst3_featureplot2.png",height = 12,width = 16)
all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 2) 
top10.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker20.csv")

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",9),9)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 3,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP.png",dpi = 300,width = 8,height = 6)
DimPlot(int.seu, group.by = "subtype",reduction = "umap", label = T,repel = T,label.size = 2,label.box = T,cols = sample(MetBrewer::met.brewer("Klimt",4),3),label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP_subtype.png",dpi = 300,width = 8,height = 6)
qs::qsave(int.seu,file = "./vst3.qs")


int.seu$NF2_status <- ifelse(int.seu$orig.ident%in%c("T3","T4","T5","T7"),"Frame-shift-del","Normal Nerve")
int.seu$NF2_status <- ifelse(int.seu$orig.ident%in%c("T2","T6"),"Nonsense",int.seu$NF2_status)
int.seu$NF2_status <- ifelse(int.seu$orig.ident=="T1","Wildtype",int.seu$NF2_status)
int.seu$NF2_status <- ifelse(int.seu$orig.ident%in%c("N1","N2"),"Normal Nerve",int.seu$NF2_status)

seurat_colors <- sample(MetBrewer::met.brewer("Klimt",9),9)
DimPlot(int.seu, split.by = "NF2_status",reduction = "umap", label = T,repel = T,label.size = 3,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.001,raster = F)
ggsave("UMAP_NF2.png",dpi = 300,width = 8,height = 6)
