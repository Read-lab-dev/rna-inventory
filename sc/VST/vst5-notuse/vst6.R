##Subcluster##
.libPaths()
setwd("~/rna/sc/VST/vst6")
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
int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 2000,verbose=F)
int.seu = ScaleData(int.seu,verbose=F)
int.seu = RunPCA(int.seu,verbose=F)
ElbowPlot(int.seu,ndims = 50)
int.seu <- harmony::RunHarmony(int.seu, group.by.vars = "orig.ident", 
                               reduction = "pca", assay.use = "RNA", 
                               reduction.save = "harmony",max_iter=20)
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu = FindNeighbors(int.seu,dims=seq(30),verbose=F,reduction = "harmony")
int.seu <- FindClusters(int.seu,resolution = 0.3)
DimPlot(int.seu,label = T)+DimPlot(int.seu,group.by = "orig.ident")

library(DoubletFinder)
library(Seurat)
library(dplyr)

sweep.res.list <- paramSweep_v3(int.seu, PCs = 1:20, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = 0.039                     # 5000细胞对应的doublets rate是3.9%
homotypic.prop <- modelHomotypic(int.seu$seurat_clusters)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(int.seu)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
int.seu <- doubletFinder_v3(int.seu, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, 
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
table(int.seu$DF.classifications_0.25_0.05_461)
DimPlot(int.seu, reduction = "umap", group.by = "DF.classifications_0.25_0.05_461")
int.seu<- subset(int.seu,subset = DF.classifications_0.25_0.05_461 =="Singlet")
qs::qsave(int.seu,"vst6.qs")
#######################
int.seu <- qs::qread("vst6.qs")
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
ggsave("vst6_featureplot.png",height = 12,width = 16)

features <-c("EPCAM","PECAM1","MME","NRXN1","ITGB8","CD3G","CD3E","CD79A","PTPRC")
FeaturePlot(int.seu,features = unlist(features),order = T,
            reduction = "umap",raster = F)&scale_colour_viridis_c(option = "B")
ggsave("vst6_featureplot2.png",height = 12,width = 16)
all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 1) 
top10.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker13.csv")

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",7),6)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 3,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP.png",dpi = 300,width = 8,height = 6)
DimPlot(int.seu, group.by = "orig.ident",reduction = "umap", label = T,repel = T,label.size = 2,label.box = T,cols = sample(MetBrewer::met.brewer("Klimt",4),4),label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP_subtype.png",dpi = 300,width = 8,height = 6)
qs::qsave(int.seu,file = "./vst6.qs")

