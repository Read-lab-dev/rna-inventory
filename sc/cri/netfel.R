rm(list = ls())
setwd("~/rna2/cri")
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(gghighlight)
library(MetBrewer)
library(tibble)
seurat_colors <- as.character(met.brewer("Klimt", 8))
colors_continuous <- as.character(rev(met.brewer("Hiroshige", 100)))
colors_scale <- as.character(rev(met.brewer("Johnson",100)))
sub_cols <- c("#CF5A79","#544799")

int.seu <- qs::qread("gsc.seu.qs")
lineage <- read.csv("./netfel/netfel.csv",header = T,na.strings = "")

int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$G1.S)),nbin = 30,ctrl = 100,name = "G1S")
int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$G2.M)),nbin = 30,ctrl = 100,name = "G2M")

figdata <- data.frame(int.seu$G1S1,int.seu$G2M1)

data<- fitdistrplus::fitdist(int.seu$G2M1, "norm") 

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

threshold <- c(quantile(gsc.module[gsc.module$state=="MES",1],probs = 0.9,names = F),
               quantile(gsc.module[gsc.module$state=="NPC",2],probs = 0.9,names = F),
               quantile(gsc.module[gsc.module$state=="AC",3],probs = 0.9,names = F),
               quantile(gsc.module[gsc.module$state=="OPC",4],probs = 0.9,names = F)
               )
names(threshold) <- c("MES","NPC","AC","OPC")

fun.minus <- function(x){as.numeric(sort(x)[3]-sort(x)[2]) } 
gsc.module$state.hybrid <-ifelse(apply(gsc.module[,1:4],1,fun.minus)>0.1,"hybrid",gsc.module$state)


fun.above <- function(x){ names(sort(x)[3]) } 
fun.above.num <-function(x){ as.numeric(sort(x)[3]) } 
gsc.module$orig.cycling <- int.seu$cycling
idx <- apply(gsc.module[gsc.module$state.hybrid=="hybrid",1:4],1,fun.above)
num.idx <- apply(gsc.module[gsc.module$state.hybrid=="hybrid",1:4],1,fun.above.num)
gsc.module[gsc.module$state.hybrid=="hybrid","state.hybrid"]<-ifelse(num.idx>threshold[idx],
                                                                     "hybrid",gsc.module[gsc.module$state.hybrid=="hybrid",]$state)

gsc.module$D <- ifelse(gsc.module$state%in%c("NPC","MES"),0,1)

yaxis <- function(x){2*(max(x[4],x[2])-max(x[3],x[1]))}

xaxis <- function(x){aa =ifelse(x[5]>0, log2(abs(x[4]-x[2])+1),
                            log2(abs(x[3]-x[1])+1))
  aa= ifelse(x[6]>0,-2*aa,2*aa)
  return(aa)
}
gsc.module$yaxis <- apply(gsc.module[,c(1:4)],1,yaxis)
gsc.module$xaxis <- apply(gsc.module[,c(1:4,9,8)], 1, xaxis)
gsc.module$orig.ident <- int.seu$orig.ident
gsc.module$cluster <- int.seu$RNA_snn_res.0.6
gsc.module$phase.own <- int.seu$phase.own
int.seu$state <- gsc.module$state
gsc.module$state.hybrid <- ifelse(gsc.module$orig.cycling=="cycling","cycling",gsc.module$state.hybrid)
int.seu$state.hybrid <- gsc.module$state.hybrid
DimPlot(int.seu,group.by = "state.hybrid",cols = seurat_colors)
# gsc.module$state <- ifelse(abs(gsc.module$xaxis)<1&abs(gsc.module$yaxis)<1,"undiff",gsc.module$state)
gsc.module$orig.ident <- factor(gsc.module$orig.ident,levels = c("DMSO","CrizS","CrizL"))

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
         theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
         theme(axis.ticks.x = element_blank(),
               legend.position = "right",
               panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
         stat_density2d(data = gsc.module[gsc.module$orig.cycling=="cycling",],
                               aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_gradientn(colors = colors_continuous)+facet_wrap(vars(orig.ident))

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
  scale_fill_gradientn(colors = colors_continuous)+facet_wrap(vars(orig.ident))


ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=state))+
  geom_point(alpha=1,size=0.5)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373","skyblue","grey"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  facet_wrap(vars(cluster))

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=state))+
  geom_point(alpha=1,size=0.5)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373","skyblue","grey"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))

ggsave("netfel_subtype.pdf",height = 4,width = 6)


ggplot(data = gsc.module %>% filter(state.hybrid!="hybrid"),aes(x=xaxis,y=yaxis,color=cycling.score))+
  geom_point(alpha=1,size=0.5)+
  scale_color_manual(values = c("#CF5A79","grey80"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))

gsc.module$MET <- log2(GetAssayData(int.seu,slot = "count")["MET",]+1)
gsc.module$eGFP <- log2(GetAssayData(int.seu,slot = "count")["eGFP",]+1)
gsc.module$CD24 <- log2(GetAssayData(int.seu,slot = "count")["CD24",]+1)
gsc.module$CD44 <- log2(GetAssayData(int.seu,slot = "count")["CD44",]+1)
gsc.module$CD133 <- log2(GetAssayData(int.seu,slot = "count")["PROM1",]+1)
gsc.module$NES <- log2(GetAssayData(int.seu,slot = "count")["NES",]+1)
gsc.module$SOX2 <- log2(GetAssayData(int.seu,slot = "count")["ALK",]+1)

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$CD44>2,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_gradientn(colors = colors_continuous)+facet_wrap(vars(orig.ident))

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("CD133+ Cell\n AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$SOX2>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F,n=200)+
  scale_fill_gradientn(colors = colors_continuous)+facet_wrap(vars(orig.ident))

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=MET))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) +
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  geom_point(aes(color=SOX2),size=0.1,alpha=1)+
  scale_color_gradient2(low = "grey",mid="white",high="purple")+facet_wrap(vars(orig.ident))

cell.prop <- as.data.frame(prop.table(table(gsc.module$celltype, gsc.module$orig.ident)))
colnames(cell.prop) <- c("cluster", "sample", "proportion")
ggplot(cell.prop, aes(sample, proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373","grey"))+
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.1, "cm"),
    legend.key.size = unit(10, "pt"),
    axis.text = element_text(size = 5, face = "bold"),
    title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_legend(title = NULL))

ggplot(data = mes.module,aes(x=xaxis,y=yaxis,color=new.state))+
  geom_point(alpha=1,size=0.5)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373","grey"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))

boxplot(MES~new.state,data = mes.module)

boxplot(MES~orig.cycling,data = mes.module)

ggplot(gsc.module,aes(orig.ident,MET))+
  geom_boxplot()+
  facet_wrap(vars(state))


int.seu$state.hybrid <- gsc.module$state.hybrid
int.seu$state <- gsc.module$state
int.seu$orig.cycling <- gsc.module$orig.cycling

qs::qsave(int.seu,file = "gsc.seu.qs")
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
stem.data <- stem.data %>% filter(orig.ident%in%c("eng_DMSO","eng_VP1","eng_VP14"))

ggboxplot(stem.data,x="orig.ident",y="YAP1",fill = "orig.ident", bxp.errorbar = T,
          color = "black",palette="npg",facet.by = "state")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        axis.title.x = element_blank())+
  stat_compare_means(label = "p.format",
                ref.group = "eng_DMSO", 
                method = "wilcox.test")
ggsave("YAP1_state.pdf",height = 4,width = 3)

ggboxplot(stem.data,x="orig.ident",y="MYC",fill = "orig.ident", bxp.errorbar = T,
          color = "black",palette="npg",facet.by = "state")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        axis.title.x = element_blank())+
  stat_compare_means(label = "p.format",
                     ref.group = "eng_DMSO", 
                     method = "wilcox.test")
ggsave("SOX2_state.pdf",height = 4,width = 3)

ggplot(gsc.module %>% filter(state=="MES"),aes(x=MES,color=orig.ident))+
  geom_density()

