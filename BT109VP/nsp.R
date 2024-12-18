##Subcluster##
rm(list=ls())
gc()
library(SoupX)
library(Seurat)
library(dplyr)
path <- grep(pattern="21047FL-89-01", dir("/home/hzg/rna/seq/ProcessedBam/BT109_VP/"),value = T)
sample_names <- c("Org", "Eng_DMSO", "Eng_VP_S", "Eng_VP_L", "BT109_DMSO", "BT109_VP")

out1 <- load10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-01/outs")
out2 <- load10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-02/outs")
out3 <- load10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-03/outs")
out4 <- load10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-04/outs")
out1 = autoEstCont(out1)
out2 = autoEstCont(out2)
out3 = autoEstCont(out3)
out4 = autoEstCont(out4)
out1 = adjustCounts(out1)
out2 = adjustCounts(out2)
out3 = adjustCounts(out3)
out4 = adjustCounts(out4)

scRNAlist <- list()
sample_names <- sample_names[5:6]
for (i in 1:2) {
  scRNAlist[[i]] <- CreateSeuratObject(counts = get(paste0("out",i)),
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)

int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)

int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.8&percent.mt<25&nCount_RNA>1000& 
                    nFeature_RNA >=250)

## DoubletFinder
library(DoubletFinder)
library(Seurat)
library(dplyr)

sweep.res.list <- paramSweep_v3(int.seu, PCs = 1:20, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = 0.039                     # 5000细胞对应的doublets rate是3.9%
homotypic.prop <- modelHomotypic(int.seu$seurat_clusters)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(int.seu)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
int.seu <- doubletFinder_v3(int.seu, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, 
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
int.seu<- subset(int.seu,subset = DF.classifications_0.25_0.03_298 =="Singlet")

library(harmony)
int.seu <- NormalizeData(int.seu)
int.seu <- FindVariableFeatures(int.seu,nfeatures = 3000)
int.seu <- ScaleData(int.seu)
int.seu = RunPCA(int.seu,verbose=F)
int.seu <- harmony::RunHarmony(int.seu,group.by.vars="group",max_iter=30)
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu = FindNeighbors(int.seu,dims=seq(30),verbose=F,reduction = "harmony")
int.seu = FindClusters(int.seu,res=0.3,verbose=F)
FeaturePlot(int.seu,features = "GFAP",order = T)
DimPlot(int.seu,reduction = "umap",label = T)
qs::qsave(int.seu,file = "nsp.seu.qs")

all.markers.np <- FindMarkers(int.seu,ident.1 = "BT109_VP", only.pos = F, min.pct = 0.25)
all.markers <- FindMarkers(gsc.seu,ident.1 = "Eng_VP_S",ident.2 = "Eng_DMSO", only.pos = F, min.pct = 0.25)
all.markers.np$gene <- rownames(all.markers.np)
all.markers$gene <- rownames(all.markers)
inter <- inner_join(all.markers,all.markers.np,by="gene")
all.markers <- FindAllMarkers(int.seu, only.pos = T,logfc.threshold = 0.5, min.pct = 0.25)
top10.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

cell.prop <- as.data.frame(prop.table(table(int.seu$nmf_cluster,int.seu$orig.ident),margin = 2))
colnames(cell.prop)<- c("celltype","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=celltype))+
  scale_fill_manual(values = seurat_colors)+
  geom_bar(stat="identity",colour="#222222",width = 0.7)+
  theme_classic()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        panel.border = element_rect(fill = NA,color = "black",linetype = "solid"),
        title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
ggsave("ns_PROP_Celltpye.pdf",dpi = 300,width = 7,height = 4)
################NMF####
library(RcppML)
mat <- int.seu@assays$RNA@scale.data
model <- RcppML::nmf(mat, k = 50)
nmf_dim <- t(predict(model,mat))
colnames(nmf_dim)<- gsub("nmf","NMF_",colnames(nmf_dim))
int.seu[["nmf"]] <- CreateDimReducObject(embeddings = nmf_dim, key = "NMF_", assay = DefaultAssay(int.seu))
int.seu <- RunUMAP(int.seu,reduction = "nmf",dims = 1:50)
model <- RcppML::nmf(mat, k = 7)
nmf_dim <- as.data.frame(t(predict(model,mat)))
int.seu$nmf_cluster <- as.numeric(apply(nmf_dim,1,which.max))
DimPlot(int.seu,reduction = "umap",group.by = "nmf_cluster")
qs::qsave(int.seu,file = "nsp.seu.qs")
################Score########
int.seu <- qs::qread("nsp.seu.qs")
gsc.seu <- qs::qread("gsc.seu.qs")
int.seu <- merge(int.seu,gsc.seu)
library(ggplot2)
library(ggthemes)
library(gghighlight)
lineage <- read.csv("./NMF/module4.csv",header = T,na.strings = "")

int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$NPC)),nbin = 30,ctrl = 100,name = "NPC")
int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$AC)),nbin = 30,ctrl = 100,name = "AC")
int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$OPC)),nbin = 30,ctrl = 100,name = "OPC")
int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$cyc)),nbin = 30,ctrl = 100,name = "cycscore")
set.seed(20230829)

lineage.score <- pmax(int.seu$NPC1,int.seu$AC1)
lineage.score.plot <- lineage.score
lineage.class <- ifelse(int.seu$NPC1>int.seu$AC1,"NPC","AC")
# lineage.score.plot[lineage.score.plot<0] <- runif(length(lineage.score.plot[lineage.score.plot<0]),0,0.15)
lineage.data <- data.frame(lineage.score,lineage.score.plot,lineage.class)
lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"] <- -lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"]
lineage.data$lineage.class <- ifelse(int.seu$OPC1>pmax(int.seu$AC1,int.seu$NPC1),"OPC",lineage.class)
# lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"]),-0.1,0.1)
# 
# lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"]),-0.1,0.1)
lineage.data$stemness <- int.seu$OPC1-lineage.data$lineage.score
lineage.data <- data.frame(lineage.data,int.seu@meta.data)
lineage.data$cycling <- ifelse(lineage.data$cycscore1>0,"cycling","non-cycling")
seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
sub_cols <- c("#CF5A79","#544799")
lineage.data$SOX2 <- GetAssayData(int.seu)["SOX2",]
lineage.data$orig.ident<- factor(lineage.data$orig.ident,
  levels = c("Eng_DMSO", "Eng_VP_S", "Eng_VP_L", "BT109_DMSO", "BT109_VP"))

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=orig.ident))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors[c(1,3,5,7,9)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident),nrow = 1)
ggsave("ns_vs_eng.png",height = 3,width = 10)

