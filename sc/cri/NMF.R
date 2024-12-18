###NMF cluster######
rm(list = ls())
setwd("/home/hzg/rna2/cri")
int.seu <- qs::qread("int.seu.qs")
library(Seurat)
library(dplyr)
library(NMF)
int.seu <- subset(int.seu,downsample=100)
int.seu <- NormalizeData(int.seu) %>% FindVariableFeatures() %>% ScaleData(do.center = F)
vm <- as.matrix(gsc.seu@assays$RNA@scale.data)
vm[vm<=0] <- 0
vm <- vm[rowSums(vm)>0,]
res.rank <- nmf(vm, 
                 rank = 2:10,
                 seed = "random",
                 nrun=5,
                 .opt="vp10",
                 method = "lee")
qs::qsave(res.rank,file = "res.rank.qs",nthreads = 10L)
plot(res.rank)
res.rank5 <- nmf(vm, 
                rank = 8,
                seed = 123,
                nrun=5,
                .opt="vp10",
                method = "lee")
fs <- extractFeatures(res.rank5, 50L)
fs <- lapply(fs, function(x) rownames(res.rank5)[x])
fs <- do.call("rbind",fs)
rownames(fs) <- paste0("cluster", seq(nrow(fs)))
write.csv(t(fs), "pb,c_NMF_TopGenes.csv")
DT::datatable(t(fs))
s.f = 1:8
seurat.color <- sample(MetBrewer::met.brewer("Klimt",n=nlevels(Idents(int.seu))),nlevels(Idents(int.seu)))
seurat.color <- c("#469d76","#c93f55","#3c4b99","#df9ed4","#eacc62","#924099")
DimPlot(int.seu,pt.size = 0.1,cols = seurat.color)+theme(legend.position = "right",legend.text = element_text(size=6))&NoAxes()
cell.prop<-as.data.frame(prop.table(table(int.seu$state, int.seu$orig.ident)))
colnames(cell.prop)<- c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  scale_fill_manual(values = seurat.color)+
  geom_bar(stat="identity",position="fill")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1),
                   legend.position = "left")
ggsave(p1+p2,file="UMAP_NMF1.pdf",height = 4,width = 8)

qs::qsave(res.rank5,file = "res.rank5.qs",nthreads = 10L)
lineage <- read.csv("./nmf_c8.csv",header = T,na.strings = "")[1:30,]
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
  gsc.module$MES <- apply(gsc.module[,c(1,2)],1,max)
  gsc.module$NPC <- apply(gsc.module[,c(5,6)],1,max)
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
  return(int.seu)
}
int.seu <- calculate_state(int.seu)
DimPlot(int.seu,reduction = "dim2",group.by = "state",cols = seurat.color[c(2,1,3,5)])
qs::qsave(int.seu,file = "tumor.qs",nthreads = 10L)

##########Downstream#######
module <- function(int.seu){
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
gsc.module$MES <- apply(gsc.module[,c(1,2)],1,max)
gsc.module$NPC <- apply(gsc.module[,c(5,6)],1,max)
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
return(gsc.module)
}
dat <- module(int.seu)
dat$ident <- factor(dat$ident,levels = c("DMSO","cz_short","cz_long"))
dat$MET <- GetAssayData(int.seu)["SAMD5",]

FeaturePlot(int.seu,features = "DCX",order = T,reduction = "dim2")

ggplot(data = dat,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = dat[dat$MET>1,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_viridis_c(option = "E",name="MET density")+
  scale_alpha_continuous(range = c(0,1),guide=guide_none())+
  facet_wrap(vars(ident))
