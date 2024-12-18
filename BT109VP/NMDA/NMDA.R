library(Seurat)
setwd("~/rna/BT109VP/")
int.seu <- qs::qread("int.seu.NEW1.qs")
genelist <- c("GRIN2A","GRIN2B","GRIN2C","GRIN2D","GRIA1","GRIA2","GRIA3","GRIA4")
DimPlot(int.seu,label = T)+DotPlot(int.seu,features = genelist)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("BT109-ORG-NMDA.png",height = 6,width = 14)

nsp.seu <- qs::qread("../nsp.seu.qs")
DimPlot(nsp.seu,label = T)+DotPlot(nsp.seu,features = genelist)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("BT109-ns-NMDA.png",height = 6,width = 14)

tarun.seu <- qs::qread("~/rna/sc/tarun/int.seu.qs")
DimPlot(tarun.seu,label = T)+DotPlot(tarun.seu,features = genelist)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("Tarun-ns-NMDA.png",height = 6,width = 14)

gsc.seu <- qs::qread("../gsc.seu.qs")
DimPlot(gsc.seu,reduction = "dim2",label = T)+DotPlot(gsc.seu,features = genelist)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("BT109-GSC-NMDA.png",height = 6,width = 14)
