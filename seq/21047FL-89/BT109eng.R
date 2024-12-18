rm(list = ls())
setwd("~/rna/ProcessedBam/BT109_VP")
library(Seurat)
library(ggplot2)
library(dplyr)

samples <- list.files()
sample_name <- c("org", "eng_DMSO", "eng_VP1", "eng_VP14", "BT109_DMSO", "BT109_VP")
seurat_list <- list()

for (sample in samples) {
  data.path <- paste0(sample, "/outs/filtered_feature_bc_matrix")
  seurat_data <- Read10X(data.dir = data.path)
  seurat_obj <- CreateSeuratObject(
    counts = seurat_data,
    project = sample,
    min.features = 200,
    min.cells = 3
  )
  seurat_list <- append(seurat_list, seurat_obj)
}
int.seu <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = samples)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA) / log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI"), ncol = 4)

FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

int.seu <- subset(int.seu, subset = log10GenesPerUMI > 0.85 & percent.mt < 20 &
  nFeature_RNA >= 250 & nCount_RNA >= 750)
# Standard PCA
int.seu <- int.seu %>%
  Seurat::NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfLeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = F)

ElbowPlot(int.seu, ndims = 50, reduction = "pca")
pc.num <- 1:20

int.seu <- int.seu %>%
  RunUMAP(reduction = "pca", dims = pc.num, verbose = F) %>%
  FindNeighbors(reduction = "pca", k.param = 10, dims = pc.num)
DimPlot(int.seu)
library(clustree)
obj <- FindClusters(int.seu, resolution = seq(0.3, 1.0, by = 0.1))
clustree(obj)
rm(obj)
DimPlot(int.seu)

int.seu <- int.seu %>%
  FindClusters(resolution = 0.3) %>%
  identity()

all.markers <- FindAllMarkers(int.seu,
  only.pos = TRUE, min.pct = 0.25,
  logfc.threshold = 0.25, max.cells.per.ident = 300
)

top10.markers <- all.markers[-c(
  grep(c("^RP[SL]"), rownames(all.markers)),
  grep(c("^MT-"), rownames(all.markers)),
  grep(c("^ENSG"), rownames(all.markers))
), ] %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <- top10.markers %>%
    filter(cluster == i) %>%
    dplyr::select(gene) %>%
    as.data.frame()

  cluster <- paste(tmp$gene, collapse = ",")

  cellmarker <- rbind(cellmarker, cluster)
}
rownames(cellmarker) <- paste0("cluster", levels(int.seu))
write.csv(cellmarker, file = "cellmarker.csv")

new.cluster.ids <- read.csv("cellmarker.csv", header = T)$V2 ## Copy from your handmade-Celltype
Idents(int.seu) <- int.seu$seurat_clusters
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
int.seu$celltype <- factor(int.seu$celltype, levels = c(
  "GSC1", "GSC2", "VIM+RGC", "TOP2A+NPC", "SOX11+NPC", 
  "APOE+N", "ADARB2+N", "NEUROD1+N", "RELN+N",
  "NEUROD6+IPs", "A/Oligo", "FTH1+Oligo", "BCL11B+A"
))
Idents(int.seu) <- int.seu$celltype

library(MetBrewer)
seurat_colors <- as.character(met.brewer("Klimt", 13))
DimPlot(int.seu,
  reduction = "umap", label = T, repel = T, label.size = 4,
  label.box = T, label.color = "black", raster = T, raster.dpi = c(1000, 1000),
  pt.size = 2.5, cols = seurat_colors
) + NoLegend() + labs(x = "UMAP1", y = "UMAP2") +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"))
ggsave(filename = "UMAP.pdf", height = 6.5, width = 8.5)

int.seu$orig.ident <- gsub("21047FL-89-01-01", "org", int.seu$orig.ident)
int.seu$orig.ident <- gsub("21047FL-89-01-02", "eng_DMSO", int.seu$orig.ident)
int.seu$orig.ident <- gsub("21047FL-89-01-03", "eng_VP1", int.seu$orig.ident)
int.seu$orig.ident <- gsub("21047FL-89-01-04", "eng_VP14", int.seu$orig.ident)
int.seu$orig.ident <- gsub("21047FL-89-01-05", "BT109_DMSO", int.seu$orig.ident)
int.seu$orig.ident <- gsub("21047FL-89-01-06", "BT109_VP", int.seu$orig.ident)
int.seu$orig.ident <- factor(int.seu$orig.ident)
table(int.seu$orig.ident)

cell.prop <- as.data.frame(prop.table(table(Idents(int.seu), int.seu$orig.ident)))

colnames(cell.prop) <- c("cluster", "sample", "proportion")

ggplot(cell.prop, aes(sample, proportion, fill = cluster)) +
  scale_fill_manual(values = seurat_colors) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.1, "cm"),
    legend.key.size = unit(10, "pt"),
    axis.text = element_text(size = 5, face = "bold"),
    title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_legend(title = NULL))
ggsave("cell.prop.pdf", width = 4, height = 3)

saveRDS(int.seu, file = "int.seu.Rds")
