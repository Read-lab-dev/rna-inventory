rm(list = ls())
setwd("/home/hzg/rna/sc/GSE182109/")

int.seu <- qs::qread("/home/hzg/rna/sc/GSE182109/Abdelfattah2022-10x.seu.qs")
int.seu <- subset(int.seu,celltype%in%c("Tumor","Tumor-cyc"))
lineage <- read.csv("/home/hzg/rna/sc/cri/netfel/netfel.csv",header = T,na.strings = "")
calculate_state <- function(int.seu){
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
  library(ggsci)
  library(gghighlight)
  library(MetBrewer)
  library(tibble)
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
  gsc.module$D <- ifelse(gsc.module$state%in%c("NPC","MES"),0,1)
  yaxis <- function(x){2*(max(x[4],x[2])-max(x[3],x[1]))}
  
  xaxis <- function(x){aa =ifelse(x[5]>0, log2(abs(x[4]-x[2])+1),
                                  log2(abs(x[3]-x[1])+1))
  aa= ifelse(x[6]>0,-2*aa,2*aa)
  return(aa)
  }
  gsc.module$yaxis <- apply(gsc.module[,c(1:4)],1,yaxis)
  gsc.module$xaxis <- apply(gsc.module[,c(1:4,7,6)], 1, xaxis)
  gsc.module$ident <- int.seu$orig.ident
  int.seu$state <- gsc.module$state
  redc <- as.matrix(gsc.module[,8:7])
  colnames(redc) <- c("Dim2_1","Dim2_2")
  int.seu[["dim2"]] <- CreateDimReducObject(embeddings = redc,key = "Dim2_")
  DimPlot(int.seu,reduction = "dim2")
  return(int.seu)
}
int.seu <- calculate_state(int.seu)
ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=state))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
deg_volcano()
