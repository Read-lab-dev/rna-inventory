library(Seurat)
gsc.seu <-qs::qread("gsc.seu.qs")
gsc.seu$cyc <-ifelse(gsc.seu$celltype=="tCYC","cycling","non-cycling")
gsc.seu$cyc <- factor(gsc.seu$cyc,levels = c("non-cycling","cycling"))
seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
levels(gsc.seu$orig.ident) <- c("Eng_DMSO","Eng_VP_S","Eng_VP_L")
DimPlot(gsc.seu,reduction = "dim2",group.by = c("celltype"))

PTPRZ1
PTN
EGFR
MYC
SOX2
YAP1
WWTR1
NRXN1
NRXN3
NLGN3
NOTCH1
NOTCH2
NOTCH3
DLL3
ggsave("2dim-exp.png",height = 12,width = 10)
gene <- clipr::read_clip()
p1=FeaturePlot(gsc.seu,features = gene[1:7],order = T,reduction = "dim2",split.by = "orig.ident",pt.size = 0.25,ncol = 3)&
  scale_color_gradientn(colors=c("grey" ,"red"))
p2=FeaturePlot(gsc.seu,features = gene[8:14],order = T,reduction = "dim2",split.by = "orig.ident",pt.size = 0.25,ncol = 3)&
  scale_color_gradientn(colors=c("grey" ,"red"))
p1|p2
ggsave("2dim-gene-by-treat.png",height = 14,width = 16)

#
H3.1
H3C1
H3C2
H3C3
H3C4
H3C6
H3C7
H3C8
H3C10
H3C11
H3C12
H3-3A
H3-3B
gene <- clipr::read_clip()

p1=FeaturePlot(gsc.seu,features = gene[1:6],order = T,reduction = "dim2",pt.size = 0.25,split.by = "orig.ident")&
  scale_color_gradientn(colors=c("grey" ,"red"))
p2=FeaturePlot(gsc.seu,features = gene[7:12],order = T,reduction = "dim2",pt.size = 0.25,split.by = "orig.ident")&
  scale_color_gradientn(colors=c("grey" ,"red"))
p1|p2
ggsave("2dim-HIST-by-treat.png",height = 12,width = 16)

