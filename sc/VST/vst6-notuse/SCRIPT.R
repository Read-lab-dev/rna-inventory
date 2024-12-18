##Subcluster##
.libPaths()
setwd("~/rna/sc/vst2")
rm(list=ls())
gc()
library(Seurat)
library(dplyr)
packageVersion("Seurat")
path <- dir("./raw/")
sample_names <- path
scRNAlist <- list()
int.seu <- subset(int.seu,cells = umap$cellname)
for (i in 1:length(path)) {
  tmp <- Read10X(paste0("./raw/",path[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = tmp,
                                       project = sample_names[i],min.cells=3,min.features = 100)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
  }
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")+
   FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
int.seu <- subset(int.seu, subset = nFeature_RNA <6000 & nCount_RNA < 10000 &
                    percent.mt < 25)
int.seu = NormalizeData(int.seu,verbose=F)
int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=F)
int.seu = ScaleData(int.seu,verbose=F)
int.seu$technique <- ifelse(int.seu$orig.ident%in%c("VSF01","VSF02","VSF03"),"scRNA","snRNA")
# int.seu <- SCTransform(int.seu,vars.to.regress = c("nCount_RNA","percent.mt","technique"))
int.seu = RunPCA(int.seu,verbose=F)
ElbowPlot(int.seu,ndims = 50)
int.seu <- harmony::RunHarmony(int.seu, group.by.vars = "technique", 
                               reduction = "pca", assay.use = "RNA", 
                               reduction.save = "harmony",max_iter=20)
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu = FindNeighbors(int.seu,dims=seq(30),verbose=F,reduction = "harmony")
int.seu <- FindClusters(int.seu,resolution = 0.1,method = "igraph",algorithm = 4)
DimPlot(int.seu)
 aa <- FindAllMarkers(int.seu,max.cells.per.ident = 2000)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker_new.csv")
#######################
int.seu <- qs::qread("sch.seu.qs")
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",16),16)
DimPlot(int.seu,reduction = "umap", 
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

library(ggplot2)
Idents(int.seu) <- int.seu$RNA_snn_res.0.1
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",8),8)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP.png",dpi = 300,height = 6)

DimPlot(int.seu, group.by = "treat",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = sample(MetBrewer::met.brewer("Klimt",4),2),label.color = "grey100",
        pt.size = 0.1)

ggsave("UMAP_sham.png",dpi = 200)

qs::qsave(int.seu,file = "./sch.seu.qs")
