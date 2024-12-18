#2024-01-10
#engrafted organoid GFP+ cells vs. engrafted organoid GFP-, (DMSO, 24 hrs, 14 days) 
rm(list=ls())
setwd("~/lab/BT109VP/drug_response")
library(Seurat)
library(dplyr)
int.seu <- qs::qread("../int.qs")
int.seu <- subset(int.seu,subset = orig.ident%in%c("Eng_DMSO","Eng_VP_L","Eng_VP_S"))
standard10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  int.seu = NormalizeData(dat,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),verbose=verbose)
  return(int.seu)
}
int.seu <- standard10X(int.seu)
int.seu$eng <- GetAssayData(int.seu[["RNA"]],slot = "data")["eGFP",]
int.seu$eng <-ifelse(int.seu$eng==0,"org","tumor")

int.seu$cell_orig <- paste0(int.seu$orig.ident,"_",int.seu$celltype)
cellfordeg<-levels(int.seu$celltype)
Idents(int.seu) <- int.seu$cell_orig
#########Engrafted Long vs Engrafted DMSO############

DEG_RESULT <- NULL
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(int.seu, ident.1 = paste0("Eng_VP_L_",cellfordeg[i]), 
                         ident.2 = paste0("Eng_DMSO_",cellfordeg[i]), verbose = FALSE,logfc.threshold = 0.5,test.use = "MAST")
  CELLDEG$celltype <- cellfordeg[i]
  CELLDEG$symbol <- rownames(CELLDEG)
  DEG_RESULT <- rbind(DEG_RESULT,CELLDEG)
}
write.csv(DEG_RESULT[,c(1,2,6)],file = "DEG_ALL_LvsC.csv")
DEG_RESULT_Malignant  <- DEG_RESULT %>% filter(celltype=="Malignant")
DEG_RESULT_other <- DEG_RESULT %>% filter(celltype!="Malignant")
DEG_same <- dplyr::inner_join(DEG_RESULT_Malignant,DEG_RESULT_other,by="symbol")
DEG_diff <- dplyr::anti_join(DEG_RESULT_Malignant,DEG_RESULT_other,by="symbol")

write.csv(DEG_same,file = "DEG_same_LvsC.csv")
write.csv(DEG_diff,file = "DEG_diff_LvsC.csv")
#########Engrafted Short vs Engrafted DMSO############

DEG_RESULT <- NULL
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(int.seu, ident.1 = paste0("Eng_VP_S_",cellfordeg[i]), 
                         ident.2 = paste0("Eng_DMSO_",cellfordeg[i]), verbose = FALSE,logfc.threshold = 0.5,test.use = "MAST")
  CELLDEG$celltype <- cellfordeg[i]
  CELLDEG$symbol <- rownames(CELLDEG)
  DEG_RESULT <- rbind(DEG_RESULT,CELLDEG)
}
write.csv(DEG_RESULT[,c(1,2,6)],file = "DEG_ALL_SvsC.csv")

DEG_RESULT_Malignant  <- DEG_RESULT %>% filter(celltype=="Malignant")
DEG_RESULT_other <- DEG_RESULT %>% filter(celltype!="Malignant")
DEG_same <- dplyr::inner_join(DEG_RESULT_Malignant,DEG_RESULT_other,by="symbol")
DEG_diff <- dplyr::anti_join(DEG_RESULT_Malignant,DEG_RESULT_other,by="symbol")

write.csv(DEG_same,file = "DEG_same_SvsC.csv")
write.csv(DEG_diff,file = "DEG_diff_SvsC.csv")


#########Engrafted Long vs Engrafted Short############
DEG_RESULT_LvsS <- NULL
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(int.seu, ident.1 = paste0("Eng_VP_L_",cellfordeg[i]), 
                         ident.2 = paste0("Eng_VP_S_",cellfordeg[i]), verbose = FALSE,logfc.threshold = 0.5,test.use = "MAST")
  CELLDEG$celltype <- cellfordeg[i]
  CELLDEG$symbol <- rownames(CELLDEG)
  DEG_RESULT_LvsS <- rbind(DEG_RESULT_LvsS,CELLDEG)
}
write.csv(DEG_RESULT_LvsS[,c(1,2,6)],file = "DEG_ALL_LvsS.csv")

DEG_RESULT_Malignant  <- DEG_RESULT_LvsS %>% filter(celltype=="Malignant")
DEG_RESULT_other <- DEG_RESULT_LvsS %>% filter(celltype!="Malignant")

DEG_same <- dplyr::inner_join(DEG_RESULT_Malignant,DEG_RESULT_other,by="symbol")
DEG_diff <- dplyr::anti_join(DEG_RESULT_Malignant,DEG_RESULT_other,by="symbol")
write.csv(DEG_same,file = "DEG_same_LvsS.csv")
write.csv(DEG_diff,file = "DEG_diff_LvsS.csv")
