library(Seurat)
library(harmony)
library(dplyr)
setwd("~/rna/sc/LGG")
x10.seu <- qs::qread("10x_int.seu.sq")
ss2.seu <- qs::qread("ss2_int.seu.sq")

int.seu <- merge(x10.seu,ss2.seu)
int.seu <- JoinLayers(int.seu)
int.seu <- standard10X(int.seu,verbose = T)
standard10X = function(int.seu,nPCs=30,res=0.2,verbose=FALSE){
  int.seu = NormalizeData(int.seu)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 2000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu,group.by.vars=c("orig.ident","platform"),max.iter.harmony = 100)
  int.seu = RunUMAP(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction="harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
int.seu <- standard10X(int.seu)
DimPlot(int.seu)
aa <- FindAllMarkers(int.seu,logfc.threshold = 2,only.pos = T)
top10.markers <- aa %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
dat <- openxlsx::read.xlsx("./INFO.xlsx")
dat <- dat[1:54,]
metadata <- left_join(int.seu@meta.data,dat,by="orig.ident")
rownames(metadata) <- rownames(int.seu@meta.data)
int.seu@meta.data <- metadata
int.seu <-qs::qread("int.seu.qs")
aa <- FindAllMarkers(int.seu,logfc.threshold = 2,only.pos = T)
library(dplyr)
top10.markers <- aa[-c(grep(c("^RP[SL]"),rownames(aa)),
                                grep(c("^ENSG"),rownames(aa)),
                                grep(c("^MIR"),rownames(aa)),
                                grep(c("^LINC"),rownames(aa)),
                                grep(c("^MT-"),rownames(aa))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker_aa.csv")

library(ggplot2)
cellmarker <- read.csv("cellmarker_aa.csv")
new.cluster.ids <- cellmarker$V2
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",8),8)
DimPlot(int.seu, group.by = c("celltype","platform"),reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("umap1.png",height = 5,width = 10)
DimPlot(int.seu, group.by = c("Diagnosis","Sex"),reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors[4:5],label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("umap2.png",height = 5,width = 10)
DimPlot(int.seu, group.by = c("Grade","Location"),reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("umap3.png",height = 5,width = 10)
qs::qsave(int.seu,file = "int.seu.qs")

