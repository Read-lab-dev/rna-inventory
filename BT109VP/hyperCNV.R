library(Seurat)
int.seu <- qs::qread("../gsc.seu.qs")
cnv.score <- read.csv("cnv_score.csv",header = T)
cnv.leiden <- read.csv("cnv_leiden.csv",header = T)

int.seu$cnv <- cnv.score$cnv_score
int.seu$cnv_leiden <- ifelse(cnv.leiden$cnv_leiden%in%c(0,4,7),"HyperCNV","Normal")
FeaturePlot(int.seu,features = "cnv",reduction = "dim2",order = T)&scale_color_viridis_c(option = "C")
DimPlot(int.seu,group.by = "cnv_leiden",reduction = "dim2",
        alpha = 0.6,order = c("normal","HyperCNV"))+scale_colour_manual(values = c("#8080804D","darkred"))
table(int.seu$cnv_leiden,int.seu$orig.ident)

int.seu$orig.ident <- factor(int.seu$orig.ident,levels = c("Eng_DMSO","Eng_VP_S","Eng_VP_L"))
cell.prop <- as.data.frame(prop.table(table(int.seu$cnv_leiden,int.seu$orig.ident),margin = 2))
colnames(cell.prop)<- c("celltype","sample","proportion")


p1=DimPlot(int.seu,group.by = "cnv_leiden",reduction = "dim2",
        alpha = 0.4,order = c("normal","HyperCNV"))+NoLegend()+scale_colour_manual(values = c("#808080","darkred"))

p2=ggplot(cell.prop,aes(sample,proportion,fill=celltype))+
  scale_fill_manual(values = c("darkred","#8080804D"))+
  geom_bar(stat="identity",colour="#222222",width = 0.7)+
  theme_classic()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        panel.border = element_rect(fill = NA,color = "black",linetype = "solid"),
        title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
patchwork::wrap_plots(p1,p2,widths = 3:2)
ggsave("hyperCNV.pdf",height = 4,width = 6)
