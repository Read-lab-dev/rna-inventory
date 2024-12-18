####CMO intergration#####
setwd("/home/hzg/rna/seq/BT109sh/multi/count/")
rm(list=ls())
gc()
library(Seurat)
library(dplyr)
packageVersion("Seurat")
path <- dir("./")
sample_names <- path

scRNAlist <- list()
for (i in 1:length(path)) {
  tmp <- Read10X(paste0("./",path[i],"/count/sample_filtered_feature_bc_matrix"))
  scRNAlist[[i]] <- CreateSeuratObject(counts = tmp$`Gene Expression`,
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}

int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)

tmp <- Read10X_h5("raw_feature_bc_matrix.h5")
int.seu <- CreateSeuratObject(counts = tmp$`Gene Expression`,min.cells=3,min.features = 200)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)
FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")+
  FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.8&percent.mt<20& 
                    nFeature_RNA >=250&nCount_RNA<50000)

standard10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  int.seu = NormalizeData(dat,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
int.seu <- CellCycleScoring(int.seu,s.features = s.genes,g2m.features = g2m.genes)
SCT10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  int.seu <- SCTransform(int.seu,vars.to.regress =  c("percent.mt","S.Score","G2M.Score"),verbose = F)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
int.seu <- qs::qread("int.seu.qs")

int.seu <- subset(int.seu,orig.ident%in%c("ENG_GFP","ENG_shYAP","ORG"))

int.seu <- standard10X(int.seu, nPCs=30, res=0.4)

DimPlot(int.seu,cells.highlight = WhichCells(int.seu,expression=eGFP>0))+DimPlot(int.seu)


all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 1) 
top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  
  cluster <- paste(tmp$gene,collapse = ",")
  
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker19.csv")

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",13),13)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)

gsc.seu <- subset()
qs::qsave(int.seu,"eng.seu")
