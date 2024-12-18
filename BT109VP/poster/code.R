
int.seu <- subset(int.seu,orig.ident%in%c("Org","Eng_DMSO"))
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",8),8)
DimPlot(int.seu, group.by = "celltype",split.by = "orig.ident",reduction = "umap", label = F,repel = T,label.size = 4,label.box = F,cols = seurat_colors,label.color = "black",
        pt.size = 0.1)&NoAxes()
ggsave("compare1.png",height = 6,width = 12)

cell.prop <- as.data.frame(prop.table(table(int.seu$celltype,factor(int.seu$orig.ident)),margin = 2))
cell.prop$sample <- factor(cell.prop$sample)
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
gsc.seu <- qs::qread("gsc.seu.qs")

gsc.seu <- qs::qread("neurosphere.qs")
gsc.seu <- subset(gsc.seu,orig.ident=="BT109DMSO")
gsc.seu <- NormalizeData(gsc.seu)
gsc.seu <- FindVariableFeatures(gsc.seu,nfeatures = 3000)
gsc.seu <- ScaleData(gsc.seu)
gsc.seu = RunPCA(gsc.seu,verbose=F)
gsc.seu <- RunUMAP(gsc.seu,reduction = "pca",dims = 1:30)

DimPlot(gsc.seu)+ggtitle("Neurosphere")
