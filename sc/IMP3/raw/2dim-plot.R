rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(gghighlight)
library(MetBrewer)
library(tibble)
setwd("./raw")
x10.seu1 <- ReadMtx("IDHwtGBM.processed.10X.counts.mtx",cells = "cells1.proc.tsv",features = "genes1.tsv")
x10.seu2 <- ReadMtx("IDHwtGBM.processed.10X.counts.2.mtx",cells = "cells2.proc.tsv",features = "genes2.tsv")
x10.seu1 <- CreateSeuratObject(counts = x10.seu1,min.cells=3,min.features = 200)
x10.seu2 <- CreateSeuratObject(counts = x10.seu2,min.cells=3,min.features = 200)
x10.seu <- merge(x10.seu1,x10.seu2)
rm(x10.seu1,x10.seu2)
####################
smart.seq <- data.table::fread("IDHwtGBM.processed.SS2.logTPM.txt.gz")
smart.seq <- as.data.frame(smart.seq)
smart.seq <- tibble::column_to_rownames(smart.seq,var = "GENE")
smart.seq <- smart.seq[,sort(colnames(smart.seq))]
int.seu <- CreateSeuratObject(smart.seq)
metadata <- read.delim2("IDHwt.GBM.Metadata.SS2.txt",header = T,row.names = 1)[-1,]
tsne <- read.delim2("IDHwtGBM.tSNE.SS2.txt",row.names = 1)[-1,]
metadata <- metadata[sort(rownames(metadata)),]
tsne <- tsne[sort(rownames(tsne)),]
identical(rownames(metadata),colnames(int.seu))
int.seu@meta.data <- cbind(int.seu@meta.data,metadata)
Idents(int.seu) <- int.seu$Sample
colnames(tsne) <- c("tSNE_1","tSNE_2","celltype")
tsne$tSNE_1 <- as.numeric(tsne$tSNE_1)
tsne$tSNE_2 <- as.numeric(tsne$tSNE_2)
tsne_emb<- as.matrix(tsne[,1:2])
int.seu[["tsne"]] <- CreateDimReducObject(tsne_emb,key = "tSNE_")
DimPlot(int.seu)

hierarchy.result <- read.delim2("IDHwt.GBM.Hierarchy.SS2.txt")[-1,]
hierarchy.result$X <- as.numeric(hierarchy.result$X)
hierarchy.result$Y <- as.numeric(hierarchy.result$Y)
attach(hierarchy.result)
hierarchy.result$cellstate <- 
  ifelse(X>0&Y>0,"NPC-like",
    ifelse(X<0&Y>0,"OPC-like",
         ifelse(X<0&Y<0,"AC-like",
                "MES-like")))
hierarchy.result$patient <- 
  grep("MGH",unlist(strsplit(hierarchy.result$NAME,split = "-")),value = T)

hierarchy.result$patient <- substr(hierarchy.result$patient,1,6)

hierarchy.result <- hierarchy.result[grepl("MGH",hierarchy.result$NAME),]


seurat.col <- sample(met.brewer("Klimt",7),4)
hierarchy.result$RIOK2 <- GetAssayData(gbm.seu,slot = "data")["RIOK2",]
hierarchy.result$IMP3 <- GetAssayData(gbm.seu,slot = "data")["IGF2BP3",]
hierarchy.result$EGFR <- GetAssayData(gbm.seu,slot = "data")["EGFR",]
hierarchy.result$CDK4 <- GetAssayData(gbm.seu,slot = "data")["CDK4",]

setwd("../")
ggplot(data = hierarchy.result,aes(x=X,y=Y,color=cellstate))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
ggsave("2dim-all.pdf")

ggplot(data = hierarchy.result,aes(x=X,y=Y,color=cellstate))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+facet_wrap(vars(patient))
ggsave("2dim-facet.pdf")

ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$RIOK2>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_gradientn(colors = colors_continuous)
ggsave("2dim-RIOK2.pdf")

ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$RIOK2>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_gradientn(colors = colors_continuous)+
  facet_wrap(vars(patient))
ggsave("2dim-RIOK2_facet.pdf")

ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$IMP3>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_gradientn(colors = colors_continuous)
ggsave("2dim-IMP3.pdf")

ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$IMP3>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_gradientn(colors = colors_continuous)+
  facet_wrap(vars(patient))
ggsave("2dim-IMP3_facet.pdf")

ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$RIOK2>0&hierarchy.result$IMP3>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_gradientn(colors = colors_continuous)+facet_wrap(vars(patient))
ggsave("2dim-RIOK2&IMP3.pdf")

ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$MYC>0&hierarchy.result$IMP3>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_gradientn(colors = colors_continuous)
ggsave("2dim-IMP3&MYC.pdf")
save(gbm.seu,int.seu,hierarchy.result,colors_continuous,seurat.col,EGFR_patient,file = "Netfel.Rdata")

###########gbm.seu#################
gbm.seu<- int.seu[,hierarchy.result$NAME]
DimPlot(gbm.seu)
dim2 <- data.frame(dim2_1=as.numeric(hierarchy.result$X),
                   dim2_2=as.numeric(hierarchy.result$Y),
                   row.names = hierarchy.result$NAME)
dim2 <- as.matrix(dim2)
gbm.seu[["dim2"]] <- CreateDimReducObject(dim2,key = "dim2_")
DimPlot(gbm.seu,reduction = "dim2")
FeaturePlot(gbm.seu,features = "RIOK2",reduction = "dim2",split.by = "Sample")

gbm1 <- subset(gbm.seu,RIOK2>0.5&IGF2BP3>0.5)
FeatureScatter(gbm1,feature1 = "RIOK2",feature2 = "IGF2BP3")


EGFR_patient  <- c("MGH66","MGH104","MGH105","MGH106","MGH113","MGH115","MGH129","MGH143","MGH151","MGH152")


int.seu <- readRDS("../int.seu.Rds")

int.seu <- subset(int.seu,seurat_clusters== 7,invert=T)
seurat_colors <- as.character(met.brewer("Klimt", 8))
colors_continuous <- as.character(rev(met.brewer("Hiroshige", 100)))
colors_scale <- as.character(rev(met.brewer("Johnson",100)))
sub_cols <- c("#CF5A79","#544799")

lineage <- read.csv("netfel30.csv",header = T,na.strings = "")

{
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$G1.S)),nbin = 30,ctrl = 100,name = "G1S")
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$G2.M)),nbin = 30,ctrl = 100,name = "G2M")
  figdata <- data.frame(int.seu$G1S1,int.seu$G2M1)
  colnames(figdata) <- c("G1S","G2M")
  figdata$phase <- ifelse(figdata$G1S>0&figdata$G2M<0,"G1S",
                          ifelse(figdata$G1S>0&figdata$G2M>0,"G1G2",
                                 ifelse(figdata$G1S<0&figdata$G2M>0,"G2M","non")))
  
  ggplot(figdata,aes(G1S,G2M,color=phase))+
    geom_point(size=0.1)+
    scale_color_manual(values = c("#CF5A79","#D2C564","#406E89","grey50"))+
    theme_few()
  
  figdata$cycling <- ifelse(figdata$phase=="non","non","cycling")
  figdata$cycling.score <- apply(figdata[,c(1,2)],1,max)
  figdata$cycling.score <- ifelse(figdata$cycling.score>0,figdata$cycling.score,0)
  
  int.seu <- AddMetaData(int.seu,figdata$cycling,col.name = "cycling")

int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$MES2)),nbin = 30,ctrl = 100,name = "MES2")
                         int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$AC)),nbin = 30,ctrl = 100,name = "AC")
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$OPC)),nbin = 30,ctrl = 100,name = "OPC")
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$NPC1)),nbin = 30,ctrl = 100,name = "NPC1")
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$NPC2)),nbin = 30,ctrl = 100,name = "NPC2")
  
  gsc.module <- data.frame(int.seu$MES11,int.seu$MES21,int.seu$AC1,
                           int.seu$OPC1,int.seu$NPC11,int.seu$NPC21)
  colnames(gsc.module) <- c("MES1","MES2","AC","OPC","NPC1","NPC2")
  
  gsc.module$state <- apply(gsc.module,1,which.max)
  gsc.module$state <- gsub(1,"MES",gsc.module$state)
  gsc.module$state <- gsub(2,"MES",gsc.module$state)
  gsc.module$state <- gsub(3,"AC",gsc.module$state)
  gsc.module$state <- gsub(4,"OPC",gsc.module$state)
  gsc.module$state <- gsub(5,"NPC",gsc.module$state)
  gsc.module$state <- gsub(6,"NPC",gsc.module$state)
  gsc.module$MES <- apply(gsc.module[,1:2],1,max)
  gsc.module$NPC <- apply(gsc.module[,5:6],1,max)
  gsc.module <-gsc.module[,c(8,9,3,4,7)]
  
  threshold <- c(quantile(gsc.module[gsc.module$state=="MES",1],probs = 0.1,names = F),
                 quantile(gsc.module[gsc.module$state=="NPC",2],probs = 0.1,names = F),
                 quantile(gsc.module[gsc.module$state=="AC",3],probs = 0.1,names = F),
                 quantile(gsc.module[gsc.module$state=="OPC",4],probs = 0.1,names = F)
  )
  names(threshold) <- c("MES","NPC","AC","OPC")
  
  fun.minus <- function(x){ as.numeric(sort(x)[3]-sort(x)[2]) } 
  gsc.module$state.hybrid <-ifelse(apply(gsc.module[,1:4],1,fun.minus)>0.1,"hybrid",int.seu$cycling)
  
  fun.above <- function(x){ names(sort(x)[3]) } 
  fun.above.num <-function(x){ as.numeric(sort(x)[3]) } 
  gsc.module$orig.cycling <- int.seu$cycling
  idx <- apply(gsc.module[gsc.module$state.hybrid=="hybrid",1:4],1,fun.above)
  num.idx <- apply(gsc.module[gsc.module$state.hybrid=="hybrid",1:4],1,fun.above.num)
  gsc.module[gsc.module$state.hybrid=="hybrid","state.hybrid"]<-ifelse(num.idx>threshold[idx],
                                                                       "hybrid",gsc.module$orig.cycling)
  
  gsc.module$D <- ifelse(gsc.module$state%in%c("NPC","MES"),0,1)
  
  yaxis <- function(x){2*(max(x[4],x[2])-max(x[3],x[1]))}
  
  xaxis <- function(x){aa =ifelse(x[5]>0, log2(abs(x[4]-x[2])+1),
                                  log2(abs(x[3]-x[1])+1))
  aa= ifelse(x[6]>0,-2*aa,2*aa)
  return(aa)
  }
  gsc.module$yaxis <- apply(gsc.module[,c(1:4)],1,yaxis)
  gsc.module$xaxis <- apply(gsc.module[,c(1:4,9,8)], 1, xaxis)
  gsc.module$ident <- int.seu$orig.ident
}

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$state=="MES",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_gradientn(colors = colors_continuous)+facet_wrap(vars(ident))

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=state))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))

ggsave("netfel_subtype.pdf",height = 4,width = 6)

int.seu$state.hybrid <- gsc.module$state.hybrid
int.seu$state <- gsc.module$state
int.seu$orig.cycling <- gsc.module$orig.cycling
int.seu$state <- ifelse(int.seu$phase.own=="non",int.seu$state,"cycling")


saveRDS(int.seu,file = "../int.seu.Rds")
saveRDS(gsc.module,file = "../gsc.module.Rds")
###########################STEM########################################3
ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=stemness))+
  geom_point(alpha=1,size=0.5)+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  scale_color_gradientn(colors = colors_continuous)
###

library(ggpubr)
stem.data <- gsc.module
stem.data$orig.ident <- int.seu$orig.ident
stem.data$cluster <- int.seu$seurat_clusters
stem.data$YAP1 <- log2(as.numeric(int.seu@assays$RNA["YAP1",])+0.01)
stem.data$WWTR1 <- log2(as.numeric(int.seu@assays$RNA["WWTR1",])+0.01)
stem.data$TEAD1 <- log2(as.numeric(int.seu@assays$RNA["TEAD1",])+0.01)
stem.data$EGFR <- log2(as.numeric(int.seu@assays$RNA["EGFR",])+0.01)
stem.data$SOX2 <- log2(as.numeric(int.seu@assays$RNA["SOX2",])+0.01)
stem.data$SOX9 <- log2(as.numeric(int.seu@assays$RNA["SOX9",])+0.01)
stem.data$MYC <- log2(as.numeric(int.seu@assays$RNA["MKI67",])+0.01)
stem.data$MET <- log2(as.numeric(int.seu@assays$RNA["MET",])+0.01)

ggboxplot(stem.data,x="orig.ident",y="MET",fill = "orig.ident", bxp.errorbar = T,
          color = "black",palette="npg",facet.by = "state")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        axis.title.x = element_blank())+
  stat_compare_means(label = "p.format",
                     ref.group = "DMSO", 
                     method = "wilcox.test")
ggsave("MET_state.pdf",height = 8,width = 6)

ggboxplot(stem.data,x="orig.ident",y="MYC",fill = "orig.ident", bxp.errorbar = T,
          color = "black",palette="npg",facet.by = "state")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        axis.title.x = element_blank())+
  stat_compare_means(label = "p.format",
                     ref.group = "DMSO", 
                     method = "wilcox.test")
ggsave("SOX2_state.pdf",height = 4,width = 3)

cell.prop <- as.data.frame(prop.table(table(int.seu$state, int.seu$seurat_clusters)))
colnames(cell.prop) <- c("cluster", "sample", "proportion")
ggplot(cell.prop, aes(sample, proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = sample(MetBrewer::met.brewer("Klimt",5),5))+
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.1, "cm"),
    legend.key.size = unit(10, "pt"),
    axis.text = element_text(size = 5, face = "bold"),
    title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_legend(title = NULL))
