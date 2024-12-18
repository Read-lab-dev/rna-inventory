##########Filbin#########
rm(list = ls())
library(Seurat)
exp <- data.table::fread("GSE102130_K27Mproject.RSEM.vh20170621.txt.gz")
exp <- as.data.frame(exp)
rownames(exp)<- exp$Gene
exp <- exp[,-1]
dmg.seu <- CreateSeuratObject(counts = exp,min.cells=3,min.features = 200)
dmg.seu$orig.ident <- stringr::str_split(colnames(dmg.seu),'-',simplify = T)[,1]
dmg.seu$orig.ident <- stringr::str_split(dmg.seu$orig.ident,'_',simplify = T)[,1]
dmg.seu$id <- ifelse(dmg.seu$orig.ident%in%c("MUV1","MUV10","MUV5","BCH1126","BCH836","BCH869"),"K27M",ifelse(
  dmg.seu$orig.ident=="Oligo","Oligo","GBM"))
table(dmg.seu$id)
library(harmony)
standard10X = function(dat,nPCs=30,res=0.8,verbose=FALSE){
  int.seu = NormalizeData(dat,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu,group.by.vars="orig.ident")
  int.seu = RunUMAP(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
dmg.seu <- standard10X(dmg.seu)
DimPlot(dmg.seu,group.by = "orig.ident",label = T)+DimPlot(dmg.seu,label = T)+DimPlot(dmg.seu,group.by ="id")
aa <- FindAllMarkers(dmg.seu,logfc.threshold = 2)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# dmg.seu <- subset(dmg.seu,seurat_clusters%in%c("5","9","8"),invert=T)
cellmarker <- NULL
for (i in levels(dmg.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(dmg.seu))
write.csv(cellmarker,file = "cellmarker.csv")

FeaturePlot(dmg.seu,features = c("APOE","GFAP","PDGFRA","MKI67","TP53","OPCML","DLG2"),order = T)

library(ggplot2)
Idents(dmg.seu) <- dmg.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(dmg.seu)
dmg.seu <- RenameIdents(dmg.seu, new.cluster.ids)
dmg.seu$celltype <- Idents(dmg.seu)
new.levels <- sort(unique(new.cluster.ids))
dmg.seu$celltype <- factor(dmg.seu$celltype,levels = new.levels)
Idents(dmg.seu) <- dmg.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",12),12)

DimPlot(dmg.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)
qs::qsave(dmg.seu,file = "filbin.qs")

###################Intergrate###########
gsc.seu <- qs::qread("../gsc.seu.new.qs")
gsc.seu@assays$SCT <- NULL
DefaultAssay(gsc.seu)<- "RNA"
gsc.seu <- merge(gsc.seu,dmg.seu)
gsc.seu <- standard10X(gsc.seu)
DimPlot(gsc.seu,group.by ="celltype",label = T)
aa <- FindAllMarkers(dmg.seu,logfc.threshold = 2)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

score <- data.frame(gsc.seu$our.NPC1,gsc.seu$our.AC1,gsc.seu$our.OPC1,gsc.seu$our.cyc1)
aa <-apply(score, 1, which.max)
aa <- as.factor(aa)
levels(aa) <- c("NPC","AC","OPC","Cyc")
gsc.seu$sig <- aa
Idents(gsc.seu) <- gsc.seu$sig
qs::qsave(gsc.seu,file = "filbin-tumor.seu.qs")