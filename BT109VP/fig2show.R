####Fig for Present######
setwd("~/rna/BT109VP")
rm(list=ls())
gc()
library(Seurat)
library(ggplot2)
library(dplyr)
int.seu <- qs::qread("int.seu.qs")
DimPlot(int.seu,label = T,label.box = T)
int.seu$celltype <- Idents(int.seu)
Idents(int.seu) <- int.seu$seurat_clusters
int.seu@assays$SCT <- NULL
int.seu$batch <-ifelse(int.seu$orig.ident=="Eng_VP_L","batch1","batch2")
standard10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  library(harmony)
  int.seu = NormalizeData(dat,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunHarmony(int.seu, group.by.vars = "batch",max_iter=20)
  int.seu = RunUMAP(int.seu,reduction = "harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,reduction = "harmony",dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
downstream =function(dat,nPCs=30,res=1.0,verbose=FALSE){
int.seu = RunUMAP(int.seu,reduction = "harmony",dims=seq(nPCs),verbose=verbose)
int.seu = FindNeighbors(int.seu,reduction = "harmony",dims=seq(nPCs),verbose=verbose)
int.seu = FindClusters(int.seu,res=res,verbose=verbose)
return(int.seu)
}
int.seu <- standard10X(int.seu,nPCs = 15,res = 0.6)
DimPlot(int.seu,label = T,label.box = T)

all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 1) 
top10.markers <- all.markers%>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
write.csv(cellmarker,file = "cellmarker_19.csv")
int.seu <- AddModuleScore(int.seu,features = list(list1=new.cluster.ids),nbin = 30,name = "AC")
VlnPlot(int.seu,features = "AC1")

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",12),12)

DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)+
  theme_void()+
  theme(legend.position = "none",title = element_blank())
  
ggsave("umap_int.png")
qs::qsave(int.seu,"int.seu.qs")


int.seu = RunUMAP(int.seu,reduction = "harmony",dims=seq(15),verbose=T)
int.seu = FindNeighbors(int.seu,reduction = "harmony",dims=seq(15),verbose=T)
library(reticulate)
Sys.setenv(reticulate_python='/home/hzg/miniconda3/envs/r-reticulate/bin/python')
reticulate::py_module_available("leidenalg")
int.seu <- FindClusters(int.seu,resolution = 0.6,algorithm = 3,method = "igraph")

int.seu = FindClusters(int.seu,res=0.6,verbose=T,method="igraph",algorithm = 4)

DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,label.color = "grey100",
        pt.size = 0.1)

############Circ############
outside <- data.frame(table(int.seu$celltype))
outside$fraction <- outside$Freq / sum(outside$Freq)
outside$ymax <- cumsum(outside$fraction)
outside$ymin <- c(0, head(outside$ymax, n = -1))
outside$labelPosition <- (outside$ymax + outside$ymin) / 2
outside$label <- paste0(outside$Var1, "\n", paste0(round(outside$fraction,3)*100,"%"))

inside <- data.frame(table(int.seu$orig.ident))
inside$fraction <- inside$Freq / sum(inside$Freq)
inside$ymax <- cumsum(inside$fraction)
inside$ymin <- c(0, head(inside$ymax, n = -1))
inside$labelPosition <- (inside$ymax + inside$ymin) / 2
inside$label <- paste0(inside$Var1, "\n", paste0(round(inside$fraction,3)*100,"%"))

p1 <-ggplot(outside, aes(
  ymax = ymax, 
  ymin = ymin,
  xmax = 5, 
  xmin = 4.5,
  fill = Var1
)) +
  geom_rect() +
  geom_rect(data = inside,aes(
    ymax = ymax, 
    ymin = ymin,
    xmax = 4.5, 
    xmin = 4,
    fill = Var1)) +
  geom_text(x = 5.5,
            aes(y = labelPosition, label = label, color = Var1),
            size = 4.5) + 
  geom_text(data = inside,x = 3.5,
            aes(y = labelPosition, label = label, color = Var1),
            size = 4.5) + 
  scale_fill_manual(values = seurat_colors) +
  scale_color_manual(values = seurat_colors) +
  coord_polar(theta = "y") +
  xlim(c(-1, 5)) +
  theme_void() +
  theme(legend.position = "none")
p2 <- DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors[1:8],label.color = "grey100",
              pt.size = 0.1)+
  theme_void()+
  theme(legend.position = "none",title = element_blank())
p1+p2

ggsave(p1,file="./circ_gsc.png",height = 16,width = 16)
ggsave(p2,file="./umap_gsc.png",height = 8,width = 8,dpi = 200)

###########FeaturePlot##################
features <- list(
  Neuron = c("OPCML","NEUROD6","DLG2","RBFOX3","STMN2"),
  IP = c("EOMES","SOX11","FOXP2"),
  GSC = c("eGFP","EGFR","PDGFRA","OLIG1","OLIG2"),
  "RG/AC" = c("VIM","CLU","SOX2","SLC1A3","GLI3"))
int.seu$celltype <- factor(int.seu$celltype,levels = new.cluster.ids)
Idents(int.seu) <- int.seu$celltype
DotPlot(int.seu, features = features)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text  = element_text(size = 6))+
  scale_color_viridis_c(option = "E")
ggsave("dotplot_new.pdf")
DimPlot(int.seu, cells.highlight= H3K27Mcells)&theme_void()+theme(legend.position = "none")
DimPlot(int.seu, cells.highlight= gfpcells)&theme_void()+theme(legend.position = "none")

sub <- read.csv("sub.csv")
module$old <- module$Martina
module <- full_join(module,sub)
module$cured <- ifelse(is.na(module$new),module$cured,module$new)


gsc.seu<-AddModuleScore(gsc.seu,features =list(na.omit(module$cured)),nbin = 30,ctrl = 100,name="GSTEM")

FeaturePlot(gsc.seu,features = "GSTEM1",reduction = "dim2")+scale_color_viridis(option = "d")

VlnPlot(gsc.seu,features = "GSTEM1")

module <- read.csv("./GSTEM/DGC.csv")

gsc.seu<-AddModuleScore(gsc.seu,features =list(na.omit(module$cured)),nbin = 30,ctrl = 100,name="DGC")
