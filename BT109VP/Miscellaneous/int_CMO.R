###intergrated#####
library(Seurat)
int.seu <-qs::qread("int.seu.new.qs")
cmo.seu <- qs::qread("~/rna/seq/ProcessedBam/BT109_CMO/int.seu.qs")

DimPlot(int.seu,group.by = "celltype",label = T,label.box = T)
int.seu$group <- ifelse(int.seu$orig.ident%in%c("ENG_GFP","ENG_shYAP","ORG"),"Y","N")

DimPlot(int.seu,label = T,label.box = T)+DimPlot(int.seu,group.by = "orig.ident",label = T,label.box = T)
all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 2) 
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
write.csv(cellmarker,file = "cellmarker_int.csv")

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",13),13)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)
ggsave("umap_int.png")
qs::qsave(int.seu,file = "int_cmo.qs")
