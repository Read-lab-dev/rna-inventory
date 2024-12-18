rm(list = ls())
gc()
CMO.seu <- readRDS("/home/hzg/rna/seq/BT85CMO/int.seu.Rds")
gsc.seu <- Read10X_h5(filename = "sample_filtered_feature_bc_matrix.h5")
int.seu <- merge(int.seu,gsc.seu)

int.seu$seurat_clusters <- "Unengrafted"
DimPlot(int.seu,group.by = "seurat_clusters")+ggtitle("Unengrafted Neurosphere")
ggsave("/home/hzg/rna/BT109VP/poster/BT85ns_UAMP.png",height = 6,width = 7)

ns109 <- qs::qread("/home/hzg/rna/BT109VP/nsp.seu.raw.qs")
ns109$group <- "neurosphere"
DimPlot(ns109,group.by = "group")+ggtitle("Unengrafted Neurosphere")
ggsave("/home/hzg/rna/BT109VP/poster/BT109ns_UAMP.png",height = 6,width = 7)

eng_109 <- qs::qread("/home/hzg/rna/BT109VP/int.seu.NEW1.qs")
eng_109 <- subset(eng_109,orig.ident%in%c("Org","Eng_DMSO"))
FeaturePlot(eng_109,features = "eGFP",order = T)

cell.prop <- as.data.frame(prop.table(table(eng_109$celltype,factor(eng_109$orig.ident)),margin = 2))

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
