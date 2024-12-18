rm(list=ls())
library(Seurat)
int.seu <- qs::qread("int.qs")
gsc.seu <- qs::qread("gsc.seu.qs")

all.markers.int <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)
all.markers.gsc <- FindAllMarkers(gsc.seu, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25,test.use = "MAST")

top10.markers.int <- all.markers.int[-c(grep(c("^RP[SL]"),rownames(all.markers.int)),
                                grep(c("^MT-"),rownames(all.markers.int))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top10.markers.gsc <- all.markers.gsc[-c(grep(c("^RP[SL]"),rownames(all.markers.gsc)),
                                grep(c("^MT-"),rownames(all.markers.gsc))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

marker <- setdiff(top10.markers.gsc$gene,top10.markers.int)
marker <- top10.markers.gsc[top10.markers.gsc$gene%in%marker,]

FeaturePlot(gsc.seu,features = c("NTS", "OPCML", "STMN2", "CNTNAP2", "THSD7A", "NRG3", "CCSER1", "DLGAP1", "DLG2", "NAV3"))
FeaturePlot(int.seu,features = top10.markers.gsc$gene[31:40])

Idents(gsc.seu) <- paste0(gsc.seu$orig.ident,"_",Idents(gsc.seu))

AA <- FindMarkers(gsc.seu,ident.1 = "Neu-like",test.use = "MAST")
join <- inner_join(AA,resdata)
