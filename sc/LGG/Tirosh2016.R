##########Tirosh2016#########
rm(list = ls())
library(Seurat)
exp <- data.table::fread("./mis/GSE70630_OG_processed_data_v2.txt.gz")
exp <- as.data.frame(exp)
rownames(exp)<- gsub("'","",exp$V1)
exp <- exp[,-1]
exp <- exp[which(apply(exp,1,sum)!=0),]
colSums((2^exp[,1:5]-1)*10) ###Check if it is TPM


Tirosh2016.seu <- CreateSeuratObject(counts = na.omit(exp))
Tirosh2016.seu$orig.ident <- stringr::str_split(colnames(Tirosh2016.seu),'-',simplify = T)[,1]
Tirosh2016.seu$orig.ident <- stringr::str_split(Tirosh2016.seu$orig.ident,'_',simplify = T)[,1]
table(Tirosh2016.seu$orig.ident)
Tirosh2016.seu$orig.ident[-grep("MGH",Tirosh2016.seu$orig.ident)] <- paste0("MGH",Tirosh2016.seu$orig.ident[-grep("MGH",Tirosh2016.seu$orig.ident)])

Tirosh2016.seu[["percent.mt"]] <- PercentageFeatureSet(Tirosh2016.seu, pattern = "^MT")
VlnPlot(Tirosh2016.seu, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
library(harmony)
standard_SS2 = function(int.seu,nPCs=30,res=0.8,verbose=FALSE){
  LayerData(int.seu,layer = "data") = LayerData(int.seu,layer = "counts")
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu,group.by.vars="orig.ident")
  int.seu = RunUMAP(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
Tirosh2016.seu <- standard_SS2(Tirosh2016.seu,verbose = T)
DimPlot(Tirosh2016.seu,group.by = "orig.ident",label = T)+DimPlot(Tirosh2016.seu,label = T)
aa <- FindAllMarkers(Tirosh2016.seu,logfc.threshold = 2,only.pos = T)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Tirosh2016.seu <- subset(Tirosh2016.seu,seurat_clusters%in%c("5","9","8"),invert=T)
cellmarker <- NULL
for (i in levels(Tirosh2016.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(Tirosh2016.seu))
write.csv(cellmarker,file = "cellmarker.csv")

FeaturePlot(Tirosh2016.seu,features = c("APOE","GFAP","PDGFRA","MKI67","TP53","OPCML","DLG2"),order = T)

library(ggplot2)
Idents(Tirosh2016.seu) <- Tirosh2016.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(Tirosh2016.seu)
Tirosh2016.seu <- RenameIdents(Tirosh2016.seu, new.cluster.ids)
Tirosh2016.seu$celltype <- Idents(Tirosh2016.seu)
new.levels <- sort(unique(new.cluster.ids))
Tirosh2016.seu$celltype <- factor(Tirosh2016.seu$celltype,levels = new.levels)
Idents(Tirosh2016.seu) <- Tirosh2016.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",12),12)

DimPlot(Tirosh2016.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)
qs::qsave(Tirosh2016.seu,file = "Tirosh2016.seu.qs")
Tirosh2016.seu <- qs::qread("Tirosh2016.seu.qs")
