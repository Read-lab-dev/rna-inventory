rm(list = ls())
library(Seurat)
library(dplyr)
setwd("/home/hzg/rna/BT109VP/velocyto")
eng_DMSO <-qs::qread("../gsc.seu.qs")
table(int.seu$orig.ident)
eng_DMSO <- subset(int.seu,orig.ident=="Eng_DMSO")
eng_DMSO <- subset(int.seu,eGFP>0)
eng_DMSO <- subset(eng_DMSO,subset = celltype%in%c("cGSC","GSC_AS","GSC_OPC"))

standard10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  int.seu = NormalizeData(dat,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 2000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}

eng_DMSO <- standard10X(int.seu, nPCs=20, res=0.8)

setwd("/home/hzg/rna/BT109VP/velocyto")
eng_DMSO <-qs::qread("../int.seu.qs")
eng_DMSO <- subset(eng_DMSO,orig.ident=="Eng_DMSO")
aa <- paste0("21047FL-89-01-02:",gsub("-1","x",Cells(eng_DMSO)))

aa <- gsub("-1","x",Cells(eng_DMSO))
aa <- gsub("Eng_DMSO_","21047FL-89-01-02:",aa)
aa <- gsub("Eng_VP_S_","21047FL-89-01-03:",aa)
aa <- gsub("Eng_VP_L_","21047FL-89-01-04:",aa)
emb <- Embeddings(eng_DMSO, reduction = "dim2")
rownames(emb) <- aa
setwd("~/rna/BT109VP/velocyto/figures_gsc_dim2")
write.csv(aa, file = "cellID_obs.csv", row.names = F)
write.csv(emb, file = "cell_embeddings.csv")
bb <- data.frame(CellID=aa,x=Idents(eng_DMSO))
write.csv(bb, file = "clusters.csv",row.names = F)

DimPlot(eng_DMSO,group.by = "celltype",label = T)

