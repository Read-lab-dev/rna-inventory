rm(list = ls())
gc()
library(Seurat)

int.seu <- qs::qread("./vst.qs")
int.seu <- DietSeurat(int.seu)
library(AnnoProbe)
geneInfor <- annoGene(rownames(int.seu), "SYMBOL", "human")
geneInfor <- geneInfor[with(geneInfor, order(chr, start)), c(1, 4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[, 1]), ]
geneInfor <- geneInfor[!geneInfor$chr %in% c("chrY", "chrX"), ] # nolint: object_name_linter. # nolint
int.seu <- subset(int.seu,downsample=5000)
int.seu <- subset(int.seu,features = geneInfor$SYMBOL)
DefaultAssay(int.seu) <- "RNA"
library(SeuratDisk)
metadata <- int.seu@meta.data
write.csv(int.seu@meta.data, file = "metadata.csv")
write.csv(geneInfor, file = "gene.csv",row.names = F)

SaveH5Seurat(int.seu, filename = "int.seu.h5Seurat")
Convert("int.seu.h5Seurat", dest = "h5ad",overwrite = T)

####################
#####DownStream#####
####################
library(MetBrewer)
library(patchwork)
library(ggplot2)
library(dplyr)

seurat_colors <- as.character(sample(met.brewer("Klimt", 15),15))
colors_continuous <- as.character(rev(met.brewer("Hiroshige", 100)))
colors_scale <- as.character(rev(met.brewer("Johnson",100)))
sub_cols <- c("#D2C564","#CF5A79")

cnv.umap <- read.csv("cnv_umap.csv",header = T)[,-1]
cnv.score <- read.csv("cnv_score.csv",header = T)
cnv.leiden <- read.csv("cnv_leiden.csv",header = T)
int.seu <- qs::qread("../int.seu.qs")
int.seu$cnv.score <- cnv.score$cnv_score
int.seu$cnv.leiden <- as.character(cnv.leiden$cnv_leiden)
DimPlot(int.seu,group.by = "cnv.leiden")
colnames(cnv.umap) <- paste0("cnv_",1:2)
rownames(cnv.umap) <- colnames(int.seu)
cnv.umap <- as.matrix(cnv.umap)
int.seu[["cnv"]] <- CreateDimReducObject(embeddings = cnv.umap,key = "cnv_",assay = DefaultAssay(int.seu))

int.seu$tumor <- ifelse(GetAssayData(int.seu)["MET",]>0,"tumor","normal")

p1 <-DimPlot(int.seu,reduction = "cnv",cols = seurat_colors,
             label = T,label.box = T,label.size = 3,
             pt.size = 0.5,group.by = "cnv.leiden")+theme(legend.position = "right",legend.key.size = unit(5,"pt"),legend.text = element_text(size = 5))&NoAxes()&NoLegend()
p2 <-FeaturePlot(int.seu,reduction = "cnv",features = "cnv.score")+scale_color_gradientn(colors = colors_scale)&NoAxes()
p3 <-DimPlot(int.seu,reduction = "cnv",cols = sub_cols,pt.size = 0.1,group.by = "tumor")+
  theme(legend.position = "top")+labs(title = element_blank())&NoAxes()
p4 <-FeaturePlot(int.seu,reduction = "cnv",features = "eGFP",pt.size = 0.1)&NoAxes()
p1+p2+p3+p4

FeatureScatter(int.seu,feature1="eGFP",feature2 = "cnv.score")

ggsave(filename = "cnv_umap.jpg",p1+p2+p3+p4,height = 6.5,width = 8,dpi = 300)

saveRDS(int.seu,"../int.seu.Rds")

FeaturePlot(int.seu,reduction = "cnv",features = "MET")+
  FeaturePlot(int.seu,reduction = "cnv",features = "cnv.score")

FeaturePlot(int.seu,features = "cnv.score")
DimPlot(int.seu,group.by = "tumor")


FeaturePlot(int.seu,reduction = "cnv",features = "eGFP")

plot.data=int.seu@meta.data
ggplot(plot.data,aes(seurat_clusters,cnv.score/min(plot.data$cnv.score)))+
  geom_boxplot()
