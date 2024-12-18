rm(list = ls())
setwd("./netfel")
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(gghighlight)
library(MetBrewer)
library(tibble)
int.seu <- readRDS("../int.seu.Rds")

int.seu <- subset(int.seu,celltype== "Neoplasm")
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
  int.seu <- AddMetaData(int.seu,figdata$cycling.score,col.name = "cycling.score")
  int.seu <- AddMetaData(int.seu,figdata$phase,col.name = "phase.own")
  
  int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$MES1)),nbin = 30,ctrl = 100,name = "MES1")
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
  
  yaxis <- function(x){1*(max(x[4],x[2])-max(x[3],x[1]))}
  
  xaxis <- function(x){aa =ifelse(x[5]>0, log2(abs(x[4]-x[2])+1),
                                  log2(abs(x[3]-x[1])+1))
  aa= ifelse(x[6]>0,-2*aa,2*aa)
  return(aa)
  }
  gsc.module$yaxis <- apply(gsc.module[,c(1:4)],1,yaxis)
  gsc.module$xaxis <- apply(gsc.module[,c(1:4,9,8)], 1, xaxis)
  gsc.module$ident <- int.seu$orig.ident
}
redc <- as.matrix(gsc.module[,10:9])
colnames(redc) <- c("Dim2_1","Dim2_2")
int.seu[["dim2"]] <- CreateDimReducObject(embeddings = redc,key = "Dim2_")
DimPlot(int.seu,reduction = "dim2")
int.seu@meta.data <- cbind(int.seu@meta.data,gsc.module[,1:11])
qs::qsave(int.seu,file = "richard.tumor.qs")
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

###########################3
gsc.module$RIOK2 <-GetAssayData(int.seu)["RIOK2",]
gsc.module$IMP3 <-GetAssayData(int.seu)["IGF2BP3",]
gsc.module$MYC <-GetAssayData(int.seu)["MYC",]

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$RIOK2>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_viridis_c(option = "D",name="RIOK2 density")+scale_alpha_continuous(range = c(0,1),guide=guide_none())+
ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$IMP3>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_viridis_c(option = "D",name="IMP3 density")+scale_alpha_continuous(range = c(0,1),guide=guide_none())+
ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$MYC>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_viridis_c(option = "D",name="MYC density")+scale_alpha_continuous(range = c(0,1),guide=guide_none())

DimPlot(int.seu,cells.highlight = WhichCells(int.seu,expression=IMP3>0),
        reduction = "dim2",sizes.highlight=0.01)
FeaturePlot(int.seu,features = "MYC",reduction = "dim2",order = T)
VlnPlot(int.seu,features = "RIOK2",group.by = "state")
