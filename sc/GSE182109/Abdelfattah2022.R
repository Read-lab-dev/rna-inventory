##########Abdelfattah2022#########
rm(list = ls())
library(Seurat)
library(dplyr)
library(reticulate)
library(anndata)
library(Seurat) 
library(Matrix)
setwd("~/rna/sc/Abdelfattah2022/")
path <- grep("GSM",dir(),value = T)
scRNAlist <- list()
for (i in 1:length(path)) {
  testdf <- Read10X(paste0("./",path[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = testdf,project = path[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id= path[i])
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
}
# Simple Merge
Abdelfattah2022.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)

Abdelfattah2022.seu <- JoinLayers(Abdelfattah2022.seu)
table(Abdelfattah2022.seu$orig.ident)
Abdelfattah2022.seu$orig.ident <- ifelse(Abdelfattah2022.seu$orig.ident=="GSM5518638","MDAG11","MDAG3")
Abdelfattah2022.seu[["percent.mt"]] <- PercentageFeatureSet(Abdelfattah2022.seu, pattern = "^MT-")
VlnPlot(Abdelfattah2022.seu, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
Abdelfattah2022.seu <- subset(Abdelfattah2022.seu, subset = percent.mt<25&nFeature_RNA >=250&nCount_RNA>=1000)
library(harmony)
standard10X = function(int.seu,nPCs=30,res=0.3,verbose=FALSE){
  int.seu = NormalizeData(int.seu)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu,group.by.vars="orig.ident")
  int.seu = RunUMAP(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
Abdelfattah2022.seu <- standard10X(Abdelfattah2022.seu,verbose = T)
DimPlot(Abdelfattah2022.seu,group.by = "orig.ident",label = T)+DimPlot(Abdelfattah2022.seu,label = T)
FeaturePlot(Abdelfattah2022.seu,features = "STK17A",order = T)
ggsave("STK17A.png",height = 4,width = 5)
aa <- FindAllMarkers(Abdelfattah2022.seu,logfc.threshold = 2,only.pos = T)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
qs::qsave(Abdelfattah2022.seu,file = "Abdelfattah2022-10x.seu.qs")

Abdelfattah2022.seu <- qs::qread("Abdelfattah2022-10x.seu.qs")
DimPlot(Abdelfattah2022.seu)

seurat_colors <- as.character(sample(met.brewer("Klimt", 14),14))
colors_continuous <- as.character(rev(met.brewer("Hiroshige", 100)))
colors_scale <- as.character(rev(met.brewer("Johnson",100)))
sub_cols <- c("#D2C564","#CF5A79")

cnv.umap <- read.csv("cnv_umap.csv",header = T)[,-1]
cnv.score <- read.csv("cnv_score.csv",header = T)
cnv.leiden <- read.csv("cnv_leiden.csv",header = T)
Abdelfattah2022.seu$cnv.score <- cnv.score$cnv_score
Abdelfattah2022.seu$cnv.leiden <- as.character(cnv.leiden$cnv_leiden)
DimPlot(Abdelfattah2022.seu,group.by = "cnv.leiden")
colnames(cnv.umap) <- paste0("cnv_",1:2)
rownames(cnv.umap) <- colnames(Abdelfattah2022.seu)
cnv.umap <- as.matrix(cnv.umap)
Abdelfattah2022.seu[["cnv"]] <- CreateDimReducObject(embeddings = cnv.umap,key = "cnv_",assay = DefaultAssay(Abdelfattah2022.seu))


p1 <-DimPlot(Abdelfattah2022.seu,reduction = "cnv",cols = seurat_colors,
             label = T,label.box = T,label.size = 3,
             pt.size = 0.5,group.by = "cnv.leiden")+theme(legend.position = "right",legend.key.size = unit(5,"pt"),legend.text = element_text(size = 5))&NoAxes()&NoLegend()
p2 <-FeaturePlot(Abdelfattah2022.seu,reduction = "cnv",features = "cnv.score")+scale_color_gradientn(colors = colors_scale)&NoAxes()
p3 <-DimPlot(Abdelfattah2022.seu,reduction = "cnv",cols = sub_cols,pt.size = 0.1,group.by = "tumor")+
  theme(legend.position = "top")+labs(title = element_blank())&NoAxes()
p4 <-FeaturePlot(Abdelfattah2022.seu,reduction = "cnv",features = "eGFP",pt.size = 0.1)&NoAxes()
p1+p2+p3+p4

#######################
Abdelfattah2022.seu <- qs::qread("Abdelfattah2022-10x.seu.qs")

Abdelfattah2022.seu <- subset(Abdelfattah2022.seu,celltype=="Tumor")
lineage.filbin <- read.csv("./neftel.csv",header = T,na.strings = "")
lineage.filbin <-lineage.filbin[1:30,]
for (i in colnames(lineage.filbin)) {
  Abdelfattah2022.seu <- AddModuleScore(Abdelfattah2022.seu,features = list(lineage.filbin[,i]),
                            nbin = 30,ctrl = 100,name = paste0("neftel.",i))
}
lineage.score <- pmax(Abdelfattah2022.seu@meta.data[,9:16])
lineage.score.plot <- lineage.score