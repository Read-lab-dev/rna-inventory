##Subcluster##
setwd("~/rna/BT109VP/Mangena2024")
rm(list=ls())
gc()
library(Seurat)
library(dplyr)
TPMs <- data.table::fread("GSE276838_HUMAN_MODELS_TPM.txt.gz")
TPMs <- tibble::column_to_rownames(TPMs,var = "V1")
colnames(TPMs) <- gsub("[.]","-",colnames(TPMs))
metadata <- data.table::fread("GSE276838_HUMAN_MODELS_METADATA.txt.gz")
metadata <- tibble::column_to_rownames(metadata,var = "V1")

int.seu <- CreateSeuratObject(counts = TPMs,meta.data = metadata)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
int.seu@assays$RNA$data <- int.seu@assays$RNA$counts 
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)

int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.8&percent.mt<25&nCount_RNA>1000& 
                    nFeature_RNA >=250)
int.seu <- JoinLayers(int.seu)
standard10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  # int.seu = NormalizeData(int.seu)
  int.seu = FindVariableFeatures(int.seu,nfeatures = 3000)
  int.seu = ScaleData(int.seu)
  int.seu = RunPCA(int.seu,verbose=verbose)
  library(harmony)
  int.seu = RunHarmony(int.seu, group.by.vars = "orig.ident", reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),reduction = "harmony",verbose=verbose)
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),reduction = "harmony",verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
int.seu <- standard10X(int.seu,res = 0.3)

all.markers <- FindAllMarkers(int.seu, only.pos = T, min.pct = 0.25,logfc.threshold = 1)

top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^ENSG"),rownames(all.markers)),
                                grep(c("^MIR"),rownames(all.markers)),
                                grep(c("^LINC"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker_aa.csv")
# FeaturePlot(int.seu,features = c("APOE","GFAP","PDGFRA","MKI67","TP53"),order = T)
DimPlot(int.seu,label = T,label.box = T)

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",8),8)
DimPlot(int.seu,group.by = c("Cell_type","Type","CNA","Organoid_line"),
        reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,
        cols = seurat_colors,label.color = "grey100",pt.size = 0.1)
ggsave("umap.png",width = 12,height = 9)

all.markers <- FindAllMarkers(int.seu,logfc.threshold = 1,only.pos = T)
top10.markers <- all.markers[!grepl(c("^ENSG"),rownames(all.markers)),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = rev(p_val_adj))
qs::qsave(int.seu,file = "int.seu.qs")


#######Integrate#############
all.seu <- qs::qread("int.seu.qs")
org.seu <- qs::qread("../int.seu.NEW1.qs")
org.seu <- subset(org.seu,orig.ident=="Org")
int.seu <- subset(all.seu,CNA=="Non-malignant")
Idents(int.seu) <- paste0("suva_",Idents(int.seu))
Idents(org.seu) <- paste0("our_",Idents(org.seu))
int.seu <- merge(int.seu,org.seu)
standard10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  int.seu = NormalizeData(int.seu)
  int.seu = FindVariableFeatures(int.seu,nfeatures = 3000)
  int.seu = ScaleData(int.seu)
  int.seu = RunPCA(int.seu,verbose=verbose)
  library(harmony)
  int.seu = RunHarmony(int.seu, group.by.vars = "orig.ident", reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),reduction = "harmony",verbose=verbose)
  # int.seu = FindNeighbors(int.seu,dims=seq(nPCs),reduction = "harmony",verbose=verbose)
  # int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
int.seu <-standard10X(int.seu)
DimPlot(int.seu,label = T)
ggsave("umap_int_our.png",height = 6,width = 9)

# cor.data <- AggregateExpression(int.seu,group.by = "ident",features=rownames(int.seu@assays$RNA$scale.data))
# cor.data <- as.data.frame(cor.data)
# colnames(cor.data) <- gsub("RNA.","",colnames(cor.data))
# library(corrplot)
# library(ggplot2)
# p.cor <-cor(cor.data)
# p.cor <- p.cor[8:14,1:7]
# testRes <- cor.mtest(cor.data,conf.level=0.95)
# testRes$p <- testRes$p[8:14,1:7]
# pdf(file = "corrplot.jessa.pdf", height = 5,width = 5)
# corrplot(p.cor,method = 'square',p.mat = testRes$p,
#          title = "Jessa.2022.Pons",
#          sig.level = c(0.001,0.01,0.05),pch.cex = 0.9,
#          insig = "label_sig",tl.col = "black",
#          col=rev(COL2('RdYlBu', 100)))
# dev.off()
