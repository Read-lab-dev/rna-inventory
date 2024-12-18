library(Seurat)
library(dplyr)
mye.seu <- qs::qread("mye_new.seu.qs")
rm(list = ls())
setwd("~/rna/sc/VST")
int.seu <- qs::qread("vst.qs")
mye.seu <- subset(int.seu,celltype=="Myeloid cell")
mye.seu = NormalizeData(mye.seu,verbose=T,assay = "RNA")
mye.seu = FindVariableFeatures(mye.seu,selection.method = "vst",nfeatures= 2000,verbose=F)
mye.seu = ScaleData(mye.seu,verbose=F)
mye.seu = RunPCA(mye.seu,verbose=F)
ElbowPlot(mye.seu,ndims = 50)
# mye.seu$technique <- ifelse(mye.seu$orig.ident%in%c("SCH18","SCH20","SCH21","SCH22"),"frozen","fresh")
# mye.seu$technique <- ifelse(mye.seu$seurat_clusters%in%c(3,5,7),"frozen","fresh")
mye.seu <- RunHarmony(mye.seu, group.by.vars = "orig.ident",ncores=10,plot_convergence=T,
                               reduction = "pca", assay.use = "RNA", lambda=1,
                               kmeans_init_nstart=20, kmeans_init_iter_max=100,
                               reduction.save = "harmony",max_iter=30)
mye.seu = FindNeighbors(mye.seu,dims=seq(30),verbose=F,reduction = "harmony")
mye.seu <- RunUMAP(mye.seu,reduction = "harmony",dims = 1:30)
mye.seu <- FindClusters(mye.seu,resolution = 0.5)
DimPlot(mye.seu,raster = F,label=T)
all.markers <- FindAllMarkers(mye.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 1.5) 
top10.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
cellmarker <- NULL
for (i in levels(mye.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(mye.seu))
write.csv(cellmarker,file = "cellmarker_myeloid.csv")
mye.seu <- subset(mye.seu,seurat_clusters%in%c(14,12),invert=T)
DimPlot(mye.seu,raster = F,label = T)

DotPlot(mye.seu,features = c("CD14","FCGR3A","CD68","CD86","MRC1","CCR5","TFRC","ITGAM"))

library(ggplot2)
Idents(mye.seu) <- mye.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(mye.seu)
mye.seu <- RenameIdents(mye.seu, new.cluster.ids)
mye.seu$celltype <- Idents(mye.seu)
new.levels <- sort(unique(new.cluster.ids))
mye.seu$celltype <- factor(mye.seu$celltype,levels = new.levels)
Idents(mye.seu) <- mye.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt", 13),13)

# mye.seu$technique <- ifelse(mye.seu$orig.ident%in%c("SCH18","SCH20","SCH21","SCH22"),"frozen","fresh")
# mye.seu$technique <- ifelse(mye.seu$batch=="vst3","BD",mye.seu$technique)
DimPlot(mye.seu, group.by = c("seurat_clusters","batch"),reduction = "umap", label = T,repel = T,label.size = 3,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.01,raster = F)
ggsave("UMAP_macro.png",dpi = 300,width = 14,height = 6)
qs::qsave(mye.seu,file = "./mye_new.seu.qs")
mye.seu <- qs::qread("mye.seu.qs")

seurat_colors <- c("#DF9ED4","#44927A","#3E5F90","#D16284","#98B46C","#924099","#D16258","#E5BA60","#5C4699","#3F688C","#BDBF67","#448E7C")

DimPlot(mye.seu, group.by = c("seurat_clusters","batch"),reduction = "umap", label = T,repel = T,label.size = 3,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP_macro.png",dpi = 300,width = 7,height = 6)

####MARCOPHAGE####
feature_list <- list(M1=c("CCL2","CCL3","CCL5","CXCL8","CD86","IL1B","HLA-DRA","HLA-DRB1"),
                     M2a=c("ITGAM","PTPRC","CD14","MRC1","CD163","IL1R1","TNF"),
                     M2b=c("TNF","IL1B","IL1A","IL6","IL10","CCL1","SPHK1"),
                     M2c=c("IL10","TGFB1","CCL16","CCL18","CXCL13"),
                     M2d=c("IL10","IL6","TGFB1","VEGFA","CCL18","CSF1"))
mye.seu <- AddModuleScore(mye.seu,features = feature_list)
VlnPlot(mye.seu,features = c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5"),pt.size=0)

feature_list <- list(General=c("CD68","CD14","ITGAM","MRC1"),
                     Mono=c("P2RY12","CX3CR1"),
                     M1=c("TNF","IL6","NOS2","CXCL9","CXCL10"),
                     M2=c("CD163","ARG1","IL10","TGFB1","MARCO"))
DotPlot(mye.seu,features = feature_list)

ggsave("dotplot_macro.png",dpi = 300,width = 12,height = 6)

mye.seu$dataset <- substr(colnames(mye.seu),1,4)
cell.prop <- as.data.frame(prop.table(table(mye.seu$seurat_clusters,mye.seu$orig.ident),margin = 2))
colnames(cell.prop)<- c("celltype","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=celltype))+
  scale_fill_manual(values = seurat_colors)+
  geom_bar(stat="identity",colour="#222222",width = 0.7)+
  theme_classic()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        panel.border = element_rect(fill = NA,color = "black",linetype = "solid"),
        title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
ggsave("porp_macro_dataset.png",dpi = 300,width = 10,height = 6)

mye.seu$celltype <- ifelse(mye.seu$celltype%in%c("5_Mono_TEX14","2_Mono_ELMO1"),"mono",mye.seu$celltype)

mye.seu$double <- ifelse(rownames(mye.seu@meta.data)%in%WhichCells(mye.seu,expression=MRC1>0&CD86>0),"T","F")
DimPlot(mye.seu,cells.highlight = WhichCells(mye.seu,expression=MRC1>0&CD86>0),raster = F,sizes.highlight = 0.1)
