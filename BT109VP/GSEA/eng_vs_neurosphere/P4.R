#2024-01-03
#engrafted organoid GFP+ cells vs. neurosphere tumor cells, (DMSO, 24 hrs, 14 days) 
rm(list=ls())
library(Seurat)
library(dplyr)
int.seu <- qs::qread("../int.qs")

int.seu <- subset(int.seu,eGFP>0)

int.seu$treat <- ifelse(int.seu$orig.ident%in%c("BT109_DMSO","Eng_DMSO"),"DMSO","VP-Treated")

standard10X = function(dat,nPCs=30,verbose=FALSE){
  int.seu = NormalizeData(dat,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 2000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),verbose=verbose)
  return(int.seu)
}
int.seu <- standard10X(int.seu)

DimPlot(int.seu,group.by = "orig.ident")

Idents(int.seu.vp)<- int.seu.vp$orig.ident

DEG.LvsC <- DEG.LvsC %>% arrange(avg_log2FC)
DEG.LvsS <- DEG.LvsS %>% arrange(avg_log2FC)
DEG.SvsC <- DEG.SvsC %>% arrange(avg_log2FC)
write.csv(DEG.LvsC,file = "DEG.LvsC.csv")
write.csv(DEG.LvsS,file = "DEG.LvsS.csv")
write.csv(DEG.SvsC,file = "DEG.SvsC.csv")

int.seu.VP <- subset(int.seu,subset = treat=="VP-Treated")
int.seu.VP <- standard10X(int.seu.VP)
Idents(int.seu.VP)<- int.seu.VP$orig.ident
DEG.LvsC <- DEG.EvsNeu <- FindMarkers(int.seu.VP,ident.1 = "Eng_VP_S",ident.2 = "BT109_VP",test.use = "MAST",logfc.threshold = 0)
DEG.LvsC <- DEG.LvsC %>% arrange(desc(avg_log2FC))
write.csv(DEG.LvsC,file = "./all/DEG.SvsC.CSV")

int.seu.DMSO <- subset(int.seu,subset = treat=="DMSO")
int.seu.DMSO <- standard10X(int.seu.DMSO)
Idents(int.seu.DMSO)<- int.seu.DMSO$orig.ident

DEG.EvsNeu <- FindMarkers(int.seu.DMSO,ident.1 = "Eng_DMSO",ident.2 = "BT109_DMSO",test.use = "MAST",logfc.threshold = 0)
DEG.EvsNeu <- DEG.EvsNeu %>% arrange(desc(avg_log2FC))
write.csv(DEG.EvsNeu,file = "./all/DEG.EvsNeu.CSV")
