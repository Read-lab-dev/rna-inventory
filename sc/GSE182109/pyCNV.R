rm(list = ls())
gc()
setwd("~/backup/lab/BT109VP/Rcnv")
library(Seurat)

Abdelfattah2022.seu <- qs::qread("../Abdelfattah2022.seu.NEW.qs")
Abdelfattah2022.seu <- DietSeurat(Abdelfattah2022.seu)
# 
library(AnnoProbe)
geneInfor <- annoGene(rownames(Abdelfattah2022.seu), "SYMBOL", "human")
geneInfor$chr <- factor(geneInfor$chr,c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX","chrY","chrM"))
geneInfor <- geneInfor[with(geneInfor, order(chr, start)), c(1, 4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[, 1]), ]
geneInfor <- geneInfor[!geneInfor$chr %in% c("chrY", "chrX"), ] # nolint: Abdelfattah2022.seu_name_linter. # nolint
Abdelfattah2022.seu <- subset(Abdelfattah2022.seu,features = geneInfor$SYMBOL)
DefaultAssay(Abdelfattah2022.seu) <- "RNA"
Abdelfattah2022.seu@assays$RNA$data <- NULL
Abdelfattah2022.seu@assays$RNA$scale.data <- NULL
library(SeuratDisk)
metadata <- Abdelfattah2022.seu@meta.data
write.csv(Abdelfattah2022.seu@meta.data, file = "metadata.csv")
write.csv(geneInfor, file = "gene.csv",row.names = F)
data.table::fwrite(Abdelfattah2022.seu@assays$RNA$counts,file = "counts.csv")

Abdelfattah2022.seu[["RNA3"]] <- as(object = Abdelfattah2022.seu[["RNA"]], Class = "Assay")
DefaultAssay(Abdelfattah2022.seu) <- "RNA3"
Abdelfattah2022.seu[["RNA"]] <- NULL
Abdelfattah2022.seu <- RenameAssays(Abdelfattah2022.seu, RNA3 = 'RNA')
SaveH5Seurat(Abdelfattah2022.seu, filename = "Abdelfattah2022.seu.h5Seurat")
Convert("Abdelfattah2022.seu.h5Seurat", dest = "h5ad",overwrite = T)

####################
#####DownStream#####
####################
library(MetBrewer)
library(patchwork)
library(ggplot2)
library(dplyr)

seurat_colors <- as.character(sample(met.brewer("Klimt", 14),14))
colors_continuous <- as.character(rev(met.brewer("Hiroshige", 100)))
colors_scale <- as.character(rev(met.brewer("Johnson",100)))
sub_cols <- c("#D2C564","#CF5A79")

cnv.umap <- read.csv("cnv_umap.csv",header = T)[,-1]
cnv.score <- read.csv("cnv_score.csv",header = T)
cnv.leiden <- read.csv("cnv_leiden.csv",header = T)
Abdelfattah2022.seu <- qs::qread("../Abdelfattah2022.seu.qs")
Abdelfattah2022.seu$cnv.score <- cnv.score$cnv_score
Abdelfattah2022.seu$cnv.leiden <- as.character(cnv.leiden$cnv_leiden)
DimPlot(Abdelfattah2022.seu,group.by = "cnv.leiden")
colnames(cnv.umap) <- paste0("cnv_",1:2)
rownames(cnv.umap) <- colnames(Abdelfattah2022.seu)
cnv.umap <- as.matrix(cnv.umap)
Abdelfattah2022.seu[["cnv"]] <- CreateDimReducAbdelfattah2022.seu(embeddings = cnv.umap,key = "cnv_",assay = DefaultAssay(Abdelfattah2022.seu))

Abdelfattah2022.seu$tumor <- ifelse(GetAssayData(Abdelfattah2022.seu)["MET",]>0,"tumor","normal")

p1 <-DimPlot(Abdelfattah2022.seu,reduction = "cnv",cols = seurat_colors,
             label = T,label.box = T,label.size = 3,
             pt.size = 0.5,group.by = "cnv.leiden")+theme(legend.position = "right",legend.key.size = unit(5,"pt"),legend.text = element_text(size = 5))&NoAxes()&NoLegend()
p2 <-FeaturePlot(Abdelfattah2022.seu,reduction = "cnv",features = "cnv.score")+scale_color_gradientn(colors = colors_scale)&NoAxes()
p3 <-DimPlot(Abdelfattah2022.seu,reduction = "cnv",cols = sub_cols,pt.size = 0.1,group.by = "tumor")+
  theme(legend.position = "top")+labs(title = element_blank())&NoAxes()
p4 <-FeaturePlot(Abdelfattah2022.seu,reduction = "cnv",features = "eGFP",pt.size = 0.1)&NoAxes()
p1+p2+p3+p4

FeatureScatter(Abdelfattah2022.seu,feature1="eGFP",feature2 = "cnv.score")

ggsave(filename = "cnv_umap.jpg",p1+p2+p3+p4,height = 6.5,width = 8,dpi = 300)

saveRDS(Abdelfattah2022.seu,"../Abdelfattah2022.seu.Rds")

FeaturePlot(Abdelfattah2022.seu,reduction = "cnv",features = "MET")+
  FeaturePlot(Abdelfattah2022.seu,reduction = "cnv",features = "cnv.score")

FeaturePlot(Abdelfattah2022.seu,features = "cnv.score")
DimPlot(Abdelfattah2022.seu,group.by = "tumor")


FeaturePlot(Abdelfattah2022.seu,reduction = "cnv",features = "eGFP")

plot.data=Abdelfattah2022.seu@meta.data
ggplot(plot.data,aes(seurat_clusters,cnv.score/min(plot.data$cnv.score)))+
  geom_boxplot()
