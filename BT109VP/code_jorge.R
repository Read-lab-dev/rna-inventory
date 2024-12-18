####Code for Jorge######
setwd("~/rna/BT109VP")
rm(list=ls())
gc()
library(Seurat)
library(ggplot2)
library(dplyr)

#Load engrafted tumor cell
gsc.seu <- qs::qread("gsc.seu.qs")

#Load all cell included organoid
int.seu <- qs::qread("int.seu.qs")

Idents(int.seu)<- paste0(int.seu$orig.ident,"_",Idents(int.seu))

DEG <- NULL
for (i in levels(int.seu$celltype)) {
  DEG_tem <- FindMarkers(int.seu,ident.1=paste0("Eng_VP_L_",i),ident.2 = paste0("Eng_DMSO_",i),test.use = "MAST")
  DEG_tem$celltype <- i
  DEG_tem <- tibble::rownames_to_column(DEG_tem,"gene")
  DEG <- rbind(DEG,DEG_tem)
}
write.csv(DEG,"DEG_VPLvsDMSO_all_celltype.csv")

DEG <- NULL
for (i in levels(int.seu$celltype)) {
  DEG_tem <- FindMarkers(int.seu,ident.1=paste0("Eng_VP_S_",i),ident.2 = paste0("Eng_DMSO_",i),test.use = "MAST")
  DEG_tem$celltype <- i
  DEG_tem <- tibble::rownames_to_column(DEG_tem,"gene")
  DEG <- rbind(DEG,DEG_tem)
}
write.csv(DEG,"DEG_VPSvsDMSO_all_celltype.csv")

DEG <- NULL
for (i in levels(int.seu$celltype)) {
  DEG_tem <- FindMarkers(int.seu,ident.1=paste0("Eng_DMSO_",i),ident.2 = paste0("Org_",i),test.use = "MAST")
  DEG_tem$celltype <- i
  DEG_tem <- tibble::rownames_to_column(DEG_tem,"gene")
  DEG <- rbind(DEG,DEG_tem)
}
DEG <- DEG %>% filter(p_val_adj<0.05)
write.csv(DEG,"DEG_DMSOvsOrg_all_celltype.csv")

makers <- FindAllMarkers(int.seu,test.use="MAST",logfc.threshold = 1)

top10.markers <- makers %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(gsc.seu,features = top10.markers$gene)


int.seu <- subset(int.seu,cells=int.seu@meta.data %>% filter(orig.ident=="Org",celltype=="GSC") %>% rownames(),invert=T)

tNPC
tAC
tOPC
tCYC
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(gsc.seu)
gsc.seu <- RenameIdents(gsc.seu, new.cluster.ids)
gsc.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
gsc.seu$celltype <- factor(gsc.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
qs::qsave(int.seu,file = "int.seu.qs")
qs::qsave(gsc.seu,file = "gsc.seu.qs")
