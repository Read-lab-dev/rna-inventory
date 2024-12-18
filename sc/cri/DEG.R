###Sub cluster######
rm(list = ls())
setwd("/home/hzg/rna2/cri")
int.seu <- qs::qread("gsc.seu.qs")

rm(list = ls())
library(Seurat)
library(dplyr)
##Load your previous generated data
DimPlot(int.seu,split.by = "orig.ident")
##Time to calculate DEG
Idents(int.seu)

int.seu$celltype <- Idents(int.seu)

int.seu$celltype.group <- paste(int.seu$celltype, int.seu$orig.ident, sep = "_")

unique(int.seu$celltype.group)

saveRDS(int.seu,file="Merged_Engrafted_BT85_14D.Rds")

Idents(int.seu) <- "celltype.group"

##Loop..
##I strongly recommend you install xlsx package cause the output sheet could be numerous.
##So xlsx Format is prefered, package xlsx requires Java Enviroment.
library(writexl)

cellfordeg<-levels(int.seu$celltype)

sheets <- NULL
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(int.seu, ident.1 = paste0(cellfordeg[i],"_CrizL"), 
                         ident.2 = paste0(cellfordeg[i],"_CrizS"), verbose = FALSE,
                         test.use = "MAST")
  CELLDEG$gene <- rownames(CELLDEG)
  sheets[[i]] <- CELLDEG
}
write_xlsx(sheets,"DEG_CELL_CrizLvsS.xlsx",)

FeaturePlot(int.seu,features = c("YAP1","TAZ","TEAD1"),
            reduction = "umap",split.by = "orig.ident",pt.size = 0.5)

FeaturePlot(int.seu,features = c("eGFP"),
            reduction = "umap",split.by = "orig.ident",pt.size = 0.5)

##Visulization
top10 <-  FindMarkers(int.seu, ident.1 = paste0(cellfordeg[2],"_VP"), ident.2 = paste0(cellfordeg[2],"_DMSO"), verbose = FALSE)  %>% top_n(n = 10, wt = abs(avg_log2FC)) %>% row.names()
top10

Idents(int.seu) <- "celltype.group"
DoHeatmap(int.seu,features = top10,size=3)

Idents(int.seu) <- "celltype"

VlnPlot(int.seu,features = top10,split.by = 'orig.ident',idents = "NPC")

##Extract expression data of Single Gend
Idents(int.seu) <- colnames(int.seu)

mymatrix <- log1p(AverageExpression(int.seu, verbose = FALSE)$RNA)

mymatrix2<-t(mymatrix)%>%as.data.frame()

mymatrix2[,1]<-int.seu$celltype

colnames(mymatrix2)[1] <- "celltype"

mymatrix2[,ncol(mymatrix2)+1]<-int.seu$orig.ident

colnames(mymatrix2)[ncol(mymatrix2)] <- "orig.ident"

library(ggplot2)

p3<- ggplot2::ggplot(mymatrix2,aes(x=celltype,y=VEGFA,fill=orig.ident))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "Expression of sox2")+
  scale_x_discrete(name="Celltype")+
  scale_fill_manual(values = c('DeepSkyBlue','Orange','pink'))
p3



# Idents(int.seu) <- "orig.ident"
# CELLDEG <- FindMarkers(int.seu, ident.1 ="VP" , ident.2 ="DMSO", verbose = FALSE)
# CELLDEG$gene <- rownames(CELLDEG)