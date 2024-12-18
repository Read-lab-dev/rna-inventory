library(Seurat)
library(dplyr)
library(ggplot2)
#########Neurosphere##########
int.seu <- qs::qread("nsp.seu.qs")
int.seu <- RunUMAP(int.seu,reduction = "pca",dims = 1:30)
int.seu = FindNeighbors(int.seu,dims=seq(30),verbose=F,reduction = "pca")
int.seu = FindClusters(int.seu,res=0.3,verbose=F)
Idents(int.seu) <- "orig.ident"
all.markers.np <- FindMarkers(int.seu,ident.1 = "BT109_VP", 
                              only.pos = F, min.pct = 0.25,test.use="MAST")
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",12),2)
DimPlot(int.seu, group.by = "orig.ident",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)+ggtitle("Neurospheres")
write.csv(all.markers.np,"./DEG/ns_DEG_VP.csv")
ggsave("./DEG/NS_UMAP.png",width = 7)

#########Organoid##########
org.seu <- qs::qread("gsc.seu.qs")
Idents(org.seu) <- "celltype"
all.markers <- FindAllMarkers(org.seu,only.pos = T, min.pct = 0.25,test.use="MAST")
write.csv(all.markers,"./DEG/eng_DEG_celltype.csv")
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",12),4)
DimPlot(org.seu, group.by = c("orig.ident","celltype"),reduction = "dim2", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.5)+ggtitle("Engrafted cells")
ggsave("./DEG/eng_UMAP.png",width = 12)

Idents(org.seu) <- paste0(org.seu$celltype,"_",org.seu$orig.ident)

all.markers <- NULL
for (i in c("tOPC","tCYC")) {
  tmp <- FindMarkers(org.seu,ident.1 = paste0(i,"_Eng_VP_L"),ident.2 = paste0(i,"_Eng_DMSO"),test.use="MAST")
  tmp$celltype <- i
  all.markers <- rbind(all.markers,tmp)
}
all.markers <- all.markers %>% filter(p_val_adj<=0.05)
write.csv(all.markers,"./DEG/eng_DEG_VP_L_vs_DMSO.csv")

all.markers <- NULL
for (i in unique(org.seu$celltype)) {
  tmp <- FindMarkers(org.seu,ident.1 = paste0(i,"_Eng_VP_S"),ident.2 = paste0(i,"_Eng_DMSO"),test.use="MAST")
  tmp$celltype <- i
  all.markers <- rbind(all.markers,tmp)
}
all.markers <- all.markers %>% filter(p_val_adj<=0.05)
write.csv(all.markers,"./DEG/eng_DEG_VP_S_vs_DMSO.csv")

all.markers <- NULL
for (i in c("tOPC","tCYC")) {
  tmp <- FindMarkers(org.seu,ident.1 = paste0(i,"_Eng_VP_L"),ident.2 = paste0(i,"_Eng_VP_S"),test.use="MAST")
  tmp$gene <- rownames(tmp)
  tmp$celltype <- i
  all.markers <- rbind(all.markers,tmp)
}
all.markers <- all.markers %>% filter(p_val_adj<=0.05)
write.csv(all.markers,"./DEG/eng_DEG_VP_L_vs_S.csv")

dmso.seu <- subset(org.seu,orig.ident=="Eng_DMSO")
dmso.deg <- FindAllMarkers(dmso.seu,test.use = "MAST")
# dmso.deg <- dmso.deg %>% filter(p_val_adj<=0.05,avg_log2FC>0)
write.csv(dmso.deg,"./DEG/DEG-celltype-in-DMSO.csv")
