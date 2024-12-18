rm(list = ls())
setwd("/home/hzg/rna2/cri")
library(Seurat)
library(ggplot2)
library(dplyr)
celltype_col <- c("#EC7F51","#C969A1", "#063E5C","#EC7E28","#528791","#F8A63A","#8E9D68","#122451","#E97B6E","#D5534D","#729684","#D9B150","#CC5265","#165E6F")
path <- grep(pattern="h5", dir(),value = T)
sample_names <- gsub(".h5","",path)
# Load and QC
scRNAlist <- list()
for (i in 1:length(path)) {
  testdf <- Read10X_h5(path[i])
  
  scRNAlist[[i]] <- CreateSeuratObject(counts = testdf,
                                       project = sample_names[i],min.cells=3,min.features = 200)
  
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id= sample_names[i])
}
BT85 <- Read10X_h5("BT85")
BT85 <- CreateSeuratObject(counts = BT85$`Gene Expression`,project = "BT85",min.cells = 3,min.features = 200)
scRNAlist[[4]] <- BT85

int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)

FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")+
  FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.85&percent.mt<20& 
                    nFeature_RNA >=250&nCount_RNA>=1000)
# Standard PCA
int.seu <- int.seu %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst",nfLeatures= 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50,verbose = F)

library(harmony)
int.seu$batch <- ifelse(int.seu@meta.data$orig.ident=="BT85","cmo","normal")
int.seu <- RunHarmony(int.seu, group.by.vars="orig.ident", max_iter = 25)

ElbowPlot(int.seu, ndims=50, reduction="pca") 
pc.num <- 1:20

int.seu <- int.seu %>% 
  RunTSNE(reduction = "pca", dims = pc.num, verbose = F) %>% 
  FindNeighbors(reduction = "pca",k.param = 10, dims = pc.num)
DimPlot(int.seu,reduction = "tsne",label = T,group.by = "state.hybrid",cols = seurat_colors)

int.seu <- int.seu %>% 
  RunUMAP(reduction = "harmony", dims = pc.num, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",k.param = 10, dims = pc.num)
DimPlot(int.seu,reduction = "umap",label = T)

library(clustree)
obj <- FindClusters(int.seu, resolution = seq(0.3,1.0,by=0.1))
clustree(obj)
rm(obj)

int.seu <- int.seu %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(int.seu, reduction = "umap", label = TRUE, pt.size = .1)

FeaturePlot(int.seu,features = "MET",reduction = "tsne")

all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25) 

top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
save(all.markers,top10.markers,file = "cluster.markers17.Rdata")

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  
  cluster <- paste(tmp$gene,collapse = ",")
  
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker17.csv")

Idents(int.seu) <- int.seu$RNA_snn_res.0.8
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- clipr::read_clip()
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype

DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = celltype_col,label.color = "grey100",
        pt.size = 0.1)+theme(legend.position = "right")&NoAxes()
DimPlot(int.seu, reduction = "tsne", label = F,repel = T,label.size = 4,group.by = "malignant",
        label.box = T,label.color = "black",raster = T,raster.dpi = c(1000,1000),
        pt.size = 2.5)+theme(legend.position = "right")

FeaturePlot(int.seu,features  = "CD246")

cell.prop<-as.data.frame(prop.table(table(int.seu$malignant, int.seu$orig.ident)))
colnames(cell.prop)<- c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  scale_fill_manual(values = rainbow(2,rev = T))+
  geom_bar(stat="identity",position="fill")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))

qs::qsave(int.seu,file = "tumor.qs")

int.seu <- qs::qread("int.seu.qs")
#############OWN###########
features <- list(
  Radialglia = c("SLC1A2","NES"),
  NPC = c("TBR2","MKI67"),
  OPC    = c("PDGFRA","CSPG4"),
  Oligo      = c("OLIG1","OLIG2"),
  Astrocyte  = c("SLC1A1","S100B","ALDH1L1"),
  MatureNeuron       = c("RBFOX3","MAP2","DLG4"),
  ImatureNeuron       = c("TBR1","DCX","NEUROD1")
)
#Sloan Lab
features <- list(
  Astrocyte = c("HES1","GFAP","AGT","GLUL"),
  Neuroepithelial = c("KLF4","VIM","VCAN"),
  Neuron    = c("RBFOX3","DCX","MAP2"),
  pNPC      = c("PAX6","SOX2","GLI3","TOP2A"),
  EarlyNeuroectodermProgenitors  = c("PCP4","NEFL","NEFM"),
  Fib       = c("COL1A1","PDGFRA","DCN"),
  GliomaStemCell= c("EGFR","CD44","PROM1","FUT4","NEU1","NES","eGFP","MET"))
DotPlot(int.seu,features = features)
DimPlot(int.seu,label = T,label.box = T)

int.seu <- subset(int.seu,seurat_clusters== 7,invert=T)
int.seu$MET.sub <- ifelse(GetAssayData(int.seu,slot = "data")["MET",],"MET+","MET-")
int.seu$GFP.sub <- ifelse(GetAssayData(int.seu,slot = "data")["eGFP",],"GFP+","GFP-")


celltype_col <- rep(celltype_col,2)
cell.prop<-as.data.frame(prop.table(table(int.seu$seurat_clusters,int.seu$orig.ident)))
colnames(cell.prop)<- c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values = celltype_col)+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))

