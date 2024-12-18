rm(list = ls())
gc()
library(Seurat)
library(dplyr)
library(ggplot2)
setwd("~/rna/sc/vst4/raw/")
path <- dir()
sample_names <- c("VS3","VS1","VS2")
scRNAlist <- NULL
for (i in 1:length(sample_names)) {
  tmp <- data.table::fread(paste0(path[i],"/",path[i],"_RSEC_MolsPerCell.csv"))
  tmp <- tibble::column_to_rownames(tmp,var = "Cell_Index")
  scRNAlist[[i]] <- CreateSeuratObject(counts = t(tmp),
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
}
setwd("~/rna/sc/vst4/")
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(scRNAlist[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
int.seu <- subset(int.seu, subset = percent.mt<25)
int.seu <- NormalizeData(int.seu)
int.seu <- FindVariableFeatures(int.seu,nfeatures = 2000)
int.seu <- ScaleData(int.seu)
int.seu <- RunPCA(int.seu,verbose=F)
int.seu <- RunUMAP(int.seu,reduction = "pca",dims = 1:25)
int.seu = FindNeighbors(int.seu,dims=seq(25),verbose=F,reduction = "harmony")
int.seu = FindClusters(int.seu,res=0.3,verbose=F)
DimPlot(int.seu,reduction = "umap",label = T)
all.markers <- FindAllMarkers(int.seu, only.pos = T,logfc.threshold = 1.5, min.pct = 0.25)
top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarkers.csv")

int.seu <- qs::qread("vst1/vst1.qs")
num <- length(levels(Idents(int.seu)))
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",num),num)
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

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype

DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 3,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP.png",dpi = 300,width = 8,height = 6)
DimPlot(int.seu, group.by = "orig.ident",reduction = "umap", label = T,repel = T,label.size = 2,label.box = T,cols = sample(MetBrewer::met.brewer("Klimt",4),3),label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP_subtype.png",dpi = 300,width = 8,height = 6)

cell.prop <- as.data.frame(prop.table(table(int.seu$celltype)))
colnames(cell.prop)<- c("sample","proportion")
label_1 <- paste0(substr(cell.prop$proportion*100,1,4),"%")
p4=ggplot(cell.prop,aes(sample,proportion,fill=sample))+
  scale_fill_manual(values = seurat_colors)+
  geom_bar(stat="identity")+
  geom_text(aes(label=label_1))+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
p4



qs::qsave(int.seu,file = "vst4.qs")



