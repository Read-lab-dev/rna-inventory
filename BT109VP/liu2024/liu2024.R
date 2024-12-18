rm(list = ls())
library(Seurat)
library(ggplot2)
library(ggthemes)
library(gghighlight)
library(MetBrewer)
setwd("~/rna/BT109VP/liu2024")
exp <- data.table::fread("./raw/GSE162989_cm_tpm_frozen_postfilter.csv.gz")
exp <- as.data.frame(exp)
rownames(exp)<- exp$V1
exp <- exp[,-1]
metadata <- read.csv("raw/GSE162989_metadata_frozen_postfilter.csv.gz",row.names = 2)
dhg.seu$celltype <- metadata$annotation
dhg.seu$orig.ident <- metadata$sample
dhg.seu <- CreateSeuratObject(counts = exp,min.cells=3,min.features = 200)
dhg.seu <- merge(dhg.seu,fdhg.seu)
# dhg.seu$orig.ident <- stringr::str_split(colnames(dhg.seu),'[.]',simplify = T)[,1]
library(harmony)
standard10X = function(dat,nPCs=30,res=0.8,verbose=FALSE){
  int.seu = NormalizeData(dat,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu,group.by.vars="orig.ident")
  int.seu = RunUMAP(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
dhg.seu <- standard10X(dhg.seu)
Idents(dhg.seu) <- dhg.seu$celltype
DimPlot(dhg.seu)
qs::qsave(dhg.seu,file = "liu2024_merge.qs")

##############################################
###########Our Dataset#############
gsc.seu <- qs::qread("liu2022.qs")
gsc.seu <- subset(gsc.seu,celltype%in%c("immune_cell","Oligodendrocytes",
             "other_non_malignant","non-malignant-other","normal_oligodendrocyte"),invert=T)
gsc.seu <- JoinLayers(gsc.seu)
lineage.our <- read.csv("./module4-liu.csv",header = T,na.strings = "")
lineage.filbin <- read.csv("./liu2022.csv",header = T,na.strings = "")
lineage.our <- lineage.our[1:30,]
lineage.filbin <-lineage.filbin[1:30,]
for (i in colnames(lineage.our)) {
  gsc.seu <- AddModuleScore(gsc.seu,features = list(lineage.our[,i]),
                            nbin = 30,ctrl = 100,name = paste0("our.",i))
}

for (i in colnames(lineage.filbin)) {
  gsc.seu <- AddModuleScore(gsc.seu,features = list(lineage.filbin[,i]),
                            nbin = 30,ctrl = 100,name = paste0("liu2022.",i))
}

cor.data <- data.frame(gsc.seu@meta.data[,9:20])
colnames(cor.data) <- substr(colnames(cor.data),1,nchar(colnames(cor.data))-1)
library(corrplot)
library(ggplot2)
p.cor <-cor(cor.data)
# p.cor <- p.cor[5:9,1:4]
testRes <- cor.mtest(cor.data,conf.level=0.95)
# testRes$p <- testRes$p[5:9,1:4]
pdf(file = "corrplot.our1.pdf", height = 5,width = 5)
corrplot(p.cor,method = 'square',p.mat = testRes$p,
         title = "Liu.2022",type = "lower",
         sig.level = c(0.001,0.01,0.05),pch.cex = 0.5,
         insig = "label_sig",tl.col = "black",
         col=rev(COL2('RdYlBu', 100)))
dev.off()

###############Plot############
for (i in colnames(lineage.our)) {
  gsc.seu <- AddModuleScore(gsc.seu,features = list(lineage.our[,i]),
                            nbin = 30,ctrl = 100,name = i)
}
lineage.score <- pmax(gsc.seu$NPC1,gsc.seu$AC1)
lineage.score.plot <- lineage.score
lineage.class <- ifelse(gsc.seu$NPC1>gsc.seu$AC1,"NPC","AC")
# lineage.score.plot[lineage.score.plot<0] <- runif(length(lineage.score.plot[lineage.score.plot<0]),0,0.15)
lineage.data <- data.frame(lineage.score,lineage.score.plot,lineage.class)
lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"] <- -lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"]
lineage.data$lineage.class <- ifelse(gsc.seu$OPC1>pmax(gsc.seu$AC1,gsc.seu$NPC1),"OPC",lineage.class)
# lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"]),-0.1,0.1)
# 
# lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"]),-0.1,0.1)
lineage.data$stemness <- gsc.seu$OPC1-lineage.data$lineage.score
lineage.data <- data.frame(lineage.data,gsc.seu@meta.data)
lineage.data$cycling <- ifelse(lineage.data$cyc1>0,"cycling","non-cycling")

seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","steelblue3","#406E89","#544799", "#924099")
sub_cols <- c("#CF5A79","#544799")

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=lineage.class))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors[c(1,3,5,7,9)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
ggsave("./2dim-class.pdf",height = 3,width = 7.5)

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color= orig.ident))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors)+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  facet_wrap(vars(orig.ident),ncol = 3)
ggsave("./2dim-patient.png",height = 7,width = 12)

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=celltype))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = met.brewer("Klimt",9))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  facet_wrap(vars(celltype))
ggsave("./2dim-celltype.png",height = 8,width = 10)

cell.prop<-as.data.frame(prop.table(table(lineage.data$lineage.class,lineage.data$location)))
colnames(cell.prop)<- c("cluster","sample","proportion")
p1=ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
    scale_fill_manual(values = met.brewer("Klimt",3))+
    geom_bar(stat="identity",position="fill")+
    theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                     title = element_blank(),
                     axis.text.x = element_text(angle = 45,hjust = 1))+
    guides(fill=guide_legend(title=NULL))
cell.prop<-as.data.frame(prop.table(table(lineage.data$lineage.class,lineage.data$age)))
colnames(cell.prop)<- c("cluster","sample","proportion")
p2=ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  scale_fill_manual(values = met.brewer("Klimt",3))+
  geom_bar(stat="identity",position="fill")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
p1+p2
ggsave("prop-location.pdf")

library(dplyr)
ggplot(data = lineage.data %>%arrange(desc(cycling)),aes(x=lineage.score.plot,y=stemness,color=cycling))+
  geom_point(alpha=0.5,size=0.8)+scale_colour_manual(values =c("darkred","gray60"))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
ggsave("2dim-cycling.pdf",height = 3,width = 7.5)
