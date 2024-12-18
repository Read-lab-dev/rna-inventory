# ###12-19-2023
# Identify tumor-cell adaptations to drug treatment in a neural microenvironment:
#  engrafted organoid GFP+ cells DMSO vs. engrafted organoid GFP+ cells VP 24 hr vs. engrafted organoid GFP+ cells VP 14 days
##
rm(list=ls())
library(Seurat)
library(dplyr)
int.seu <- qs::qread("int.qs")
gsc.seu <- subset(int.seu,eGFP>0)
Idents(gsc.seu) <- gsc.seu$orig.ident
# gsc.seu@reductions$nmf <- NULL
gsc.seu <- subset(gsc.seu,idents=c("Eng_DMSO", "Eng_VP_S", "Eng_VP_L"))
DimPlot(gsc.seu)
gsc.seu <- standard10X(gsc.seu)
library(NMF)
mat <- gsc.seu@assays$RNA@scale.data
mat[mat<=0] <- 0
mat <- mat[rowSums(mat)>0,]
geneskept1 <- ifelse(rowSums(mat, na.rm = T)>3*0.01*ncol(mat),TRUE,FALSE)
geneskept2 <- ifelse(rowSums(mat > 0, na.rm = T)>3*0.01*ncol(mat),TRUE,FALSE)
geneskept <- ifelse(geneskept1==TRUE&geneskept2==TRUE,TRUE,FALSE)
mat <- mat[geneskept,]

