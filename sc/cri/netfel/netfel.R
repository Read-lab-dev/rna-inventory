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
int.seu <- qs::qread("int.seu.qs")
int.seu <- subset(int.seu,seurat_clusters==9,invert=T)
seurat_colors <- as.character(met.brewer("Klimt", 8))
colors_continuous <- as.character(rev(met.brewer("Hiroshige", 100)))
colors_scale <- as.character(rev(met.brewer("Johnson",100)))
sub_cols <- c("#CF5A79","#544799")

lineage <- read.csv("/home/hzg/rna/sc/cri/netfel/netfel.csv",header = T,na.strings = "")

int.seu <- qs::qread("int.seu.qs")
int.seu <- subset(int.seu,orig.ident=="131")
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
gsc.module$state.hybrid <-ifelse(apply(gsc.module[,1:4],1,fun.minus)>0.1,"hybrid",gsc.module$state)

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
  scale_fill_viridis_c(option = "D",name="MES-like tumor")+
  scale_alpha_continuous(range = c(0,1),guide=guide_none())+
  facet_wrap(vars(ident))
 ggsave("netfel_MES_Cri.png",height = 6,width = 18)

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=state))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))

ggsave("netfel_subtype.png",height = 6,width = 6)

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=state))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  facet_wrap(vars(ident))

ggsave("netfel_subtype2.png",height = 6,width = 18)
int.seu$state.hybrid <- gsc.module$state.hybrid
int.seu$state <- gsc.module$state
int.seu$cycling <- gsc.module$orig.cycling
# int.seu$state <- ifelse(int.seu$phase.own=="non",int.seu$state,"cycling")


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

gsc.module$fusion <- int.seu$fusion
gsc.module$MET <- GetAssayData(int.seu)["MET",]
ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$fusion=="PTPRZ1-MET Fusion",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_viridis_c(option = "D",name="MES-like tumor")+
  scale_alpha_continuous(range = c(0,1),guide=guide_none())+
  facet_wrap(vars(ident))

ggplot(data = gsc.module %>% arrange(desc(fusion)),aes(x=xaxis,y=yaxis,color=fusion))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("grey","#CF5A79"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  facet_wrap(vars(ident))

ggplot(data = gsc.module %>% arrange(MET),aes(x=xaxis,y=yaxis,color=MET))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_viridis_c(option = "D")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  facet_wrap(vars(ident))

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$MET>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_viridis_c(option = "E",name="MET+ Cells")+
  scale_alpha_continuous(range = c(0,1),guide=guide_none())+
  ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$MET==0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_viridis_c(option = "D",name="MET- Cells")+
  scale_alpha_continuous(range = c(0,1),guide=guide_none())
ggsave("MET-GBM-distribution.png",height = 5,width = 14)
