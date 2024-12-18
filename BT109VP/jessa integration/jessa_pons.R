#######Jessa Pons############
gsc.seu <- qs::qread("/home/hzg/rna/sc/BT109_VP/Jessa/pons.seu.qs")
gsc.seu$cell <- colnames(gsc.seu)
metadata <- read.delim("/home/hzg/rna/sc/BT109_VP/Jessa/metadata.tsv")
metadata$cell <- paste0(metadata$ID_paper,"_",metadata$Cell_barcode)
meta <- inner_join(metadata,gsc.seu@meta.data,by="cell")
gsc.seu <- gsc.seu[,meta$cell]
identical(colnames(gsc.seu),meta$cell)
gsc.seu$Malignant_normal_consensus_Jessa2022 <- meta$Malignant_normal_consensus_Jessa2022
gsc.seu$celltype <- meta$Cell_type_consensus_Jessa2022
DimPlot(gsc.seu,group.by = "Malignant_normal_consensus_Jessa2022")
qs::qsave(gsc.seu,file = "jessa.seu.qs")
###########################
setwd("~/rna/BT109VP")
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gghighlight)
library(MetBrewer)
library(RColorBrewer)
library(tibble)
seurat_colors <- as.character(met.brewer("Klimt", 12))
jessa.seu <- qs::qread(file = "jessa.seu.qs")
DimPlot(jessa.seu,group.by = "Malignant_normal_consensus_Jessa2022")
jessa.seu <- subset(jessa.seu,Malignant_normal_consensus_Jessa2022=="Malignant")
jessa.seu <- subset(jessa.seu,cell=="P-6253_S-8498_TGTACAGAGTTTCTTC-1",invert=T)
lineage <- read.csv("./NMF/module4.csv",header = T,na.strings = "")
jessa.seu <- AddModuleScore(jessa.seu,features = list(na.omit(lineage$NPC)),nbin = 30,ctrl = 100,name = "NPC")
jessa.seu <- AddModuleScore(jessa.seu,features = list(na.omit(lineage$AC)),nbin = 30,ctrl = 100,name = "AC")
jessa.seu <- AddModuleScore(jessa.seu,features = list(na.omit(lineage$OPC)),nbin = 30,ctrl = 100,name = "OPC")
jessa.seu <- AddModuleScore(jessa.seu,features = list(na.omit(lineage$cyc)),nbin = 30,ctrl = 100,name = "cycscore")
set.seed(20230829)
lineage.score <- pmax(jessa.seu$NPC1,jessa.seu$AC1)
lineage.score.plot <- lineage.score
lineage.class <- ifelse(jessa.seu$NPC1>jessa.seu$AC1,"NPC","AC")
lineage.data <- data.frame(lineage.score,lineage.score.plot,lineage.class)
lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"] <- -lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"]
lineage.data$lineage.class <- ifelse(jessa.seu$OPC1>pmax(jessa.seu$AC1,jessa.seu$NPC1),"OPC",lineage.class)
lineage.data$stemness <- jessa.seu$OPC1-lineage.data$lineage.score
lineage.data <- data.frame(lineage.data,jessa.seu@meta.data)
lineage.data$cycling <- ifelse(lineage.data$cycscore1>0,"cycling","non-cycling")
seurat_colors <- met.brewer("Klimt",13)
lineage.data <- lineage.data[lineage.data$cell!="P-6253_S-8498_TGTACAGAGTTTCTTC-1",]
ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=orig.ident))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors)+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))
jessa.seu$cluster <- lineage.data$lineage.class

score <- data.frame(jessa.seu$NPC1,jessa.seu$AC1,jessa.seu$OPC1,jessa.seu$cycscore1)
aa <-apply(score, 1, which.max)
aa <- as.factor(aa)
levels(aa) <- c("NPC","AC","OPC","Cyc")
jessa.seu$sig <- aa
Idents(jessa.seu) <- jessa.seu$sig
qs::qsave(jessa.seu,file = "jessa-tumor.seu.qs")
