rm(list = ls())
gc()
library(Seurat)
scRNAlist <- list()
for (i in c(1:4)) {
  cellid <- paste0("vst",i)
  scRNAlist[[i]]  <- qs::qread(paste0(cellid,"/",cellid,".qs"))
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=cellid)
}
int.seu <- subset(int.seu,orig.ident%in%c("SCH18","SCH20","SCH21","SCH22"),invert=T)
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")+
  FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
int.seu = NormalizeData(int.seu,verbose=T,assay = "RNA")
int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 2000,verbose=F)
int.seu = ScaleData(int.seu,verbose=F)
int.seu = RunPCA(int.seu,verbose=F)
ElbowPlot(int.seu,ndims = 50)
int.seu$batch <- substr(colnames(int.seu),1,4)
int.seu <- harmony::RunHarmony(int.seu, group.by.vars = "orig.ident",
                               ncores=10,plot_convergence=T,
                               reduction = "pca", assay.use = "RNA", 
                               kmeans_init_nstart=20, kmeans_init_iter_max=100,
                               reduction.save = "harmony",max_iter=30)
int.seu = FindNeighbors(int.seu,dims=seq(30),verbose=F,reduction = "harmony")
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu <- FindClusters(int.seu,resolution = 0.1)
DimPlot(int.seu,group.by = c("seurat_clusters","celltype"),raster = F,label = T)

all.markers <- FindAllMarkers(int.seu, only.pos = TRUE,max.cells.per.ident = 3000,
                              logfc.threshold = 1.5) 
top10.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker17.csv")

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(met.brewer("Klimt",10),10)
seurat_colors <- c("#CC4E56","#DF9ED4","#44927A","#E5BA60","#5C4699","#98B46C","#924099","#D16258","#3e5f90","#D16284")

DimPlot(int.seu, group.by = c("celltype","batch"),reduction = "umap", label = T,repel = T,label.size = 3,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.01,raster = F)
ggsave("UMAP.png",dpi = 300,width = 14,height = 6)
qs::qsave(int.seu,file = "./vst.qs")

cell.prop <- as.data.frame(prop.table(table(int.seu$celltype,int.seu$orig.ident),margin = 2))
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
ggsave("PROP_Celltpye.pdf",dpi = 300,width = 7,height = 4)
