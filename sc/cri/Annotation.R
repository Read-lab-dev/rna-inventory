#Annotation
rm(list = ls())
setwd("~/rna2/cri")
library(Seurat)
library(dplyr)
int.seu <- qs::qread("gsc.seu.qs")
FeaturePlot(int.seu,features = "GFAP",order = T)
FeaturePlot(int.seu,features = "AQP4",order = T)
FeaturePlot(int.seu,features = "PDGFRA",order = T)
aa <-FindAllMarkers(int.seu,max.cells.per.ident = 500)
top10.markers <- aa %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
save(all.markers,top10.markers,file = "cluster.markers17.Rdata")

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  
  cluster <- paste(tmp$gene,collapse = ",")
  
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker171.csv")

Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- clipr::read_clip()
int.seu$celltype <- factor(int.seu$celltype,levels = unique(new.levels))
Idents(int.seu) <-int.seu$celltype

int.seu@meta.data$celltype <- as.character(int.seu@meta.data$celltype)
int.seu@meta.data$celltype[int.seu$cycling=="non"&int.seu$celltype=='cAC'] <- "imAC"
int.seu@meta.data$celltype[int.seu$celltype=='mAC'] <- "imAC"

int.seu@meta.data$celltype[GetAssayData(int.seu,slot = "count")["GFAP",]>0|GetAssayData(int.seu,slot = "count")["AQP4",]>0&int.seu$celltype=='imAC'] <- "mAC"

int.seu$celltype <- as.factor(int.seu$celltype)

Idents(int.seu)<- int.seu$celltype
DimPlot(int.seu,cols = met.brewer("Klimt",12),group.by = "celltype")
qs::qsave(int.seu,"gsc.seu.qs")

############Harmony#############
library(harmony)
bt.seu <- Read10X("/home/hzg/rna2/091423sc/cmo_old/outs/per_sample_outs/Short/count/sample_filtered_feature_bc_matrix")
bt.seu <- CreateSeuratObject(counts = bt.seu$`Gene Expression`,project = "BT85",min.cells=3,min.features = 200)
bt.seu$batch <-"b2"
int.seu$batch <-"b1"
bt.seu[["percent.mt"]] <- PercentageFeatureSet(bt.seu, pattern = "^MT-")
bt.seu$log10GenesPerUMI <- log10(bt.seu$nFeature_RNA)/log10(bt.seu$nCount_RNA)
VlnPlot(bt.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)

bt.seu <- subset(bt.seu, subset = log10GenesPerUMI>0.85&percent.mt<20& 
                    nFeature_RNA >=250&nCount_RNA>=1000)
int.seu <- merge(int.seu,bt.seu)
int.seu <- int.seu %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst",nfLeatures= 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50,verbose = F)
int.seu <- RunHarmony(int.seu, group.by.vars="batch", max_iter = 25)

ElbowPlot(int.seu, ndims=50, reduction="pca") 
pc.num <- 1:20

int.seu <- int.seu %>% 
  RunUMAP(reduction = "harmony", dims = pc.num, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",k.param = 10, dims = pc.num)
DimPlot(int.seu,label = T)
FeaturePlot(int.seu,features = c("MET","eGFP"))

############Harmony#############
library(harmony)
setwd("~/rna/scRNA_analysis/BT85CMO")
bt.seu <- readRDS("int.seu.Rds")
bt.seu <- bt.seu %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst",nfLeatures= 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50,verbose = F)
bt.seu <- bt.seu %>% 
  RunUMAP(reduction = "pca", dims = pc.num, verbose = F) %>% 
  FindNeighbors(reduction = "pca",k.param = 10, dims = pc.num)
DimPlot(bt.seu,label = T)


int.seu <- merge(int.seu,bt.seu)

int.seu <- int.seu %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst",nfLeatures= 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50,verbose = F)
int.seu <- RunHarmony(int.seu, group.by.vars="batch", max_iter = 25)

ElbowPlot(int.seu, ndims=50, reduction="pca") 
pc.num <- 1:20

int.seu <- int.seu %>% 
  RunUMAP(reduction = "harmony", dims = pc.num, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",k.param = 10, dims = pc.num)
DimPlot(int.seu,label = T)
FeaturePlot(int.seu,features = c("MET","eGFP"))