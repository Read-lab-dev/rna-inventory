library(ggplot2)
cell.prop<-as.data.frame(prop.table(table(Idents(int.seu), int.seu$orig.ident)))
colnames(cell.prop)<- c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="stack")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1),
                   legend.position = "NONE")

DotPlot(int.seu,features = c("MBP","TF", "PLP1", "MAG", "MOG", "CLDN11"))
DotPlot(int.seu,features = c("CD14","AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R"))
DotPlot(int.seu,features = c("CD44","PROM1", "FUT4", "L1CAM", "ITGA6", "PDGFRA","EGFR"))+DimPlot(int.seu,label = T,label.box = T)

DimPlot(int.seu,group.by = "state.hybrid")

av <-AverageExpression(int.seu,
                       group.by = "seurat_clusters",
                       assays = "RNA")
av=av[[1]]

cg=names(tail(sort(apply(av, 1, sd)),1000))

pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))


ggplot(int.seu@meta.data)

int.seu@meta.data$orig.ident <- gsub("aaaa","CrizS",int.seu@meta.data$orig.ident)
int.seu@meta.data$orig.ident <- gsub("CrizS","CrizL",int.seu@meta.data$orig.ident)
