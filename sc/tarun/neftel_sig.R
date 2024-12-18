library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(gghighlight)
library(MetBrewer)
library(tibble)
# gbm.seu <- qs::qread("/home/hzg/rna/sc/GSE182109/int.seu.qs")
# gbm.seu <- subset(gbm.seu,Assignment=="Glioma")
# gsc.seu <- subset(int.seu,orig.ident=="1215")
# gsc.seu <- merge(gsc.seu,gbm.seu)
# gsc.seu <- JoinLayers(gsc.seu)
# gsc.seu <- NormalizeData(gsc.seu)
gsc.seu <- qs::qread("./int.seu.qs")
lineage <-read.csv("/home/hzg/rna/sc/cri/netfel/netfel.csv",header = T,na.strings = "")
{
  gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$G1.S)[1:25]),nbin = 30,ctrl = 100,name = "G1S")
  gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$G2.M)),nbin = 30,ctrl = 100,name = "G2M")
  figdata <- data.frame(gsc.seu$G1S1,gsc.seu$G2M1)
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
  
  gsc.seu <- AddMetaData(gsc.seu,figdata$cycling,col.name = "cycling")
  gsc.seu <- AddMetaData(gsc.seu,figdata$cycling.score,col.name = "cycling.score")
  gsc.seu <- AddMetaData(gsc.seu,figdata$phase,col.name = "phase.own")
  
  gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$MES1)),nbin = 30,ctrl = 100,name = "MES1")
  gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$MES2)),nbin = 30,ctrl = 100,name = "MES2")
  gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$AC)),nbin = 30,ctrl = 100,name = "AC")
  gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$OPC)),nbin = 30,ctrl = 100,name = "OPC")
  gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$NPC1)),nbin = 30,ctrl = 100,name = "NPC1")
  gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$NPC2)),nbin = 30,ctrl = 100,name = "NPC2")
  
  gsc.module <- data.frame(gsc.seu$MES11,gsc.seu$MES21,gsc.seu$AC1,
                           gsc.seu$OPC1,gsc.seu$NPC11,gsc.seu$NPC21)
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
  gsc.module$state.hybrid <-ifelse(apply(gsc.module[,1:4],1,fun.minus)>0.1,"hybrid",gsc.module$state)
  
  fun.above <- function(x){ names(sort(x)[3]) } 
  fun.above.num <-function(x){ as.numeric(sort(x)[3]) } 
  gsc.module$orig.cycling <- gsc.seu$cycling
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
  gsc.module$ident <- gsc.seu$celltype
  gsc.module$cycling <- figdata$cycling
}
dim2 <- gsc.module[,c(10,9)]
identical(rownames(dim2),colnames(gsc.seu))
colnames(dim2) <- c("dim2_1","dim2_2")
dim2$dim2_1 <- as.numeric(dim2$dim2_1)
dim2$dim2_2 <- as.numeric(dim2$dim2_2)
dim2_emb<- as.matrix(dim2)
int.seu[["dim2"]] <- CreateDimReducObject(dim2_emb,key = "dim2_")
int.seu$state <- gsc.module$state
DimPlot(int.seu,reduction = "dim2",group.by = "state")

DimPlot(gsc.seu,label = T)+
ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=state))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  annotate("text", x = -Inf, y = -Inf, label = "AC-like", hjust = 0, vjust = 0,size=rel(3),color="#544799") +
  annotate("text", x = -Inf, y = Inf, label = "OPC-like", hjust = 0, vjust = 1,size=rel(3),color="#5DA373") +
  annotate("text", x = Inf, y = -Inf, label = "MES-like", hjust = 1, vjust = 0,size=rel(3),color="#CF5A79") +
  annotate("text", x = Inf, y = Inf, label = "NPC-like", hjust = 1, vjust = 1,size=rel(3),color="#D2C564")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  facet_wrap(vars(ident))+
ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  annotate("label", x = -Inf, y = -Inf, label = "AC-like", hjust = 0, vjust = 0,size=rel(4),color="#544799") +
  annotate("label", x = -Inf, y = Inf, label = "OPC-like", hjust = 0, vjust = 1,size=rel(4),color="#5DA373") +
  annotate("label", x = Inf, y = -Inf, label = "MES-like", hjust = 1, vjust = 0,size=rel(4),color="#CF5A79") +
  annotate("label", x = Inf, y = Inf, label = "NPC-like", hjust = 1, vjust = 1,size=rel(4),color="#D2C564")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$cycling=="cycling",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_viridis_c(option = "A",name="Cycling tumor")+
  scale_alpha_continuous(range = c(0,1),guide=guide_none())

ggsave("netfel_1219.png",height = 5,width = 18)


################
gsc.module$YAP1<- GetAssayData(int.seu,slot = "data")["YAP1",]
gsc.module$TAZ<- GetAssayData(int.seu,slot = "data")["WWTR1",]
gsc.module$TEAD1<- GetAssayData(int.seu,slot = "data")["TEAD1",]
gsc.module$TEAD2<- GetAssayData(int.seu,slot = "data")["TEAD2",]
gsc.module$TEAD3<- GetAssayData(int.seu,slot = "data")["TEAD3",]
gsc.module$EGFR<- GetAssayData(int.seu,slot = "data")["EGFR",]
ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  annotate("label", x = -Inf, y = -Inf, label = "AC-like", hjust = 0, vjust = 0,size=rel(4),color="#544799") +
  annotate("label", x = -Inf, y = Inf, label = "OPC-like", hjust = 0, vjust = 1,size=rel(4),color="#5DA373") +
  annotate("label", x = Inf, y = -Inf, label = "MES-like", hjust = 1, vjust = 0,size=rel(4),color="#CF5A79") +
  annotate("label", x = Inf, y = Inf, label = "NPC-like", hjust = 1, vjust = 1,size=rel(4),color="#D2C564")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$YAP1==0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_viridis_c(option = "A",name="Cycling tumor")+
  scale_alpha_continuous(range = c(0,1),guide=guide_none())

gsc.seu$state <- gsc.module$state
ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=state))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  annotate("text", x = -Inf, y = -Inf, label = "AC-like", hjust = 0, vjust = 0,size=rel(3),color="#544799") +
  annotate("text", x = -Inf, y = Inf, label = "OPC-like", hjust = 0, vjust = 1,size=rel(3),color="#5DA373") +
  annotate("text", x = Inf, y = -Inf, label = "MES-like", hjust = 1, vjust = 0,size=rel(3),color="#CF5A79") +
  annotate("text", x = Inf, y = Inf, label = "NPC-like", hjust = 1, vjust = 1,size=rel(3),color="#D2C564")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  DotPlot(gsc.seu,features = genelist)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("Tarun-2dim-NMDA.png",height = 6,width = 14)

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=state))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  annotate("text", x = -Inf, y = -Inf, label = "AC-like", hjust = 0, vjust = 0,size=rel(3),color="#544799") +
  annotate("text", x = -Inf, y = Inf, label = "OPC-like", hjust = 0, vjust = 1,size=rel(3),color="#5DA373") +
  annotate("text", x = Inf, y = -Inf, label = "MES-like", hjust = 1, vjust = 0,size=rel(3),color="#CF5A79") +
  annotate("text", x = Inf, y = Inf, label = "NPC-like", hjust = 1, vjust = 1,size=rel(3),color="#D2C564")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  facet_wrap(vars(ident))
ggsave("Tarun-2dim.png",height =6,width = 6)
