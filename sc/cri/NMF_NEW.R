lineage <- read.csv("./netfel/netfel.csv",header = T,na.strings = "")
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

DimPlot(int.seu,group.by = c("state","nmf_cluster"),reduction = "dim2")

int.seu <- CellCycleScoring(int.seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=state))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))

qs::qread(int.seu)

gene_matrix <- int.seu@reductions$nmf@feature.loadings %>% as.data.frame()
cell_matrix <- int.seu@reductions$nmf@cell.embeddings %>% as.data.frame()
top_50_genes <- apply(gene_matrix, 2, function(x) {
  names(sort(x, decreasing = TRUE)[1:100])
})
int.seu$nmf_cluster <- apply(cell_matrix, 1, function(x) {
  which.max(x)
})

write.csv(top_50_genes,file = "own6.csv")
all.markers <- FindAllMarkers(int.seu,logfc.threshold = 2,only.pos = T)
top10.markers <- all.markers %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DotPlot(int.seu,features = c("eGFP","MET","PTPRZ1"))+DimPlot(int.seu,label = T)

seurat_colors <- sample(MetBrewer::met.brewer("Klimt",7),7)
DimPlot(int.seu, group.by = "seurat_clusters",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)
library(ggplot2)
cell.prop <- as.data.frame(prop.table(table(int.seu$state,int.seu$orig.ident),margin = 2))
colnames(cell.prop)<- c("celltype","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=celltype))+
  scale_fill_manual(values =  c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_bar(stat="identity",colour="#222222",width = 0.7)+
  theme_classic()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        panel.border = element_rect(fill = NA,color = "black",linetype = "solid"),
        title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))

DEG <- NULL
Idents(int.seu) <-paste0(int.seu$orig.ident,"_",int.seu$state)
for (i in unique(int.seu$state)) {
  tmp <- FindMarkers(int.seu,ident.1 = paste0("cz_long","_",i),
                     ident.2 = paste0("DMSO","_",i),logfc.threshold=0.5)
  tmp$celltype <- i
  tmp$gene <- rownames(tmp)
  tmp <- filter(tmp,p_val_adj<0.05)
  DEG <- rbind(DEG,tmp)
}
options(scipen = 2)
Idents(int.seu) <- factor(Idents(int.seu),levels = c("DMSO_OPC","cz_short_OPC","cz_long_OPC",
                                                     "DMSO_NPC","cz_short_NPC","cz_long_NPC",
                                                     "DMSO_AC","cz_short_AC","cz_long_AC",
                                                     "DMSO_MES","cz_short_MES","cz_long_MES"))
int.seu$compare <- Idents(int.seu)
DotPlot(int.seu,features = c("NCAM2","SOX9"))
write.csv(DEG,file = "DEG_NEW.csv")
DotPlot(int.seu,features = c("MET","EGFR","PDGFRA","BRAF","KRAS","SOX9","YAP1","WWTR1","FTH1","FTL"),scale.by = "size")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_color_viridis_c(option = "E")
Idents(int.seu) <- int.seu$state
DimPlot(int.seu)
ggsave("dotplot.pdf",height = 4,width = 8)
qs::qsave(int.seu,"tumor.qs")
save(gsc.seu,res.rank5,file = "res.rank5")
qs::qsave(int.seu,"DMSO-tumor.qs")
