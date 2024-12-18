rm(list = ls())
options(future.globals.maxSize = 8000 * 1024^2)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
packageVersion("Seurat")
setwd("~/rna/visium/spatial")
vis.seu <- Load10X_Spatial(data.dir = "~/rna/visium/spatial/visiumHD/",bin.size = c(8,16)) 
DefaultAssay(vis.seu) <- "Spatial.008um"

vln.plot <- VlnPlot(vis.seu, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(vis.seu, features = "nCount_Spatial.008um") + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot

# normalize both 8um and 16um bins
DefaultAssay(vis.seu) <- "Spatial.008um"
vis.seu <- NormalizeData(vis.seu)

DefaultAssay(vis.seu) <- "Spatial.016um"
vis.seu <- NormalizeData(vis.seu)

# switch spatial resolution to 16um from 8um
DefaultAssay(vis.seu) <- "Spatial.016um"
p1 <- SpatialFeaturePlot(vis.seu, features = "Rorb") + ggtitle("Rorb expression (16um)")

# switch back to 8um
DefaultAssay(vis.seu) <- "Spatial.008um"
p2 <- SpatialFeaturePlot(vis.seu, features = "Hpca") + ggtitle("Hpca expression (8um)")

p1 | p2

# note that data is already normalized
DefaultAssay(vis.seu) <- "Spatial.008um"
vis.seu <- FindVariableFeatures(vis.seu)
vis.seu <- ScaleData(vis.seu)
# we select 50,0000 cells and create a new 'sketch' assay
vis.seu <- SketchData(
  object = vis.seu,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switch analysis to sketched cells
DefaultAssay(vis.seu) <- "sketch"

# perform clustering workflow
vis.seu <- FindVariableFeatures(vis.seu)
vis.seu <- ScaleData(vis.seu)
vis.seu <- RunPCA(vis.seu, assay = "sketch", reduction.name = "pca.sketch")
vis.seu <- FindNeighbors(vis.seu, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
vis.seu <- FindClusters(vis.seu, cluster.name = "seurat_cluster.sketched", resolution = 3)
vis.seu <- RunUMAP(vis.seu, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

vis.seu <- ProjectData(
  object = vis.seu,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(vis.seu) <- "sketch"
Idents(vis.seu) <- "seurat_cluster.sketched"
p1 <- DimPlot(vis.seu, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(vis.seu) <- "Spatial.008um"
Idents(vis.seu) <- "seurat_cluster.projected"
p2 <- DimPlot(vis.seu, reduction = "full.umap.sketch", label = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2

SpatialDimPlot(vis.seu, label = T, repel = T, label.size = 4)

Idents(vis.seu) <- "seurat_cluster.projected"
cells <- CellsByIdentities(vis.seu, idents = c(0, 4, 32, 34, 35))
p <- SpatialDimPlot(vis.seu,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p

# Crete downsampled object to make visualization either
DefaultAssay(vis.seu) <- "Spatial.008um"
Idents(vis.seu) <- "seurat_cluster.projected"
vis.seu_subset <- subset(vis.seu, cells = Cells(vis.seu[["Spatial.008um"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(vis.seu_subset) <- "Spatial.008um"
Idents(vis.seu_subset) <- "seurat_cluster.projected"
vis.seu_subset <- BuildClusterTree(vis.seu_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(vis.seu_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

vis.seu_subset <- ScaleData(vis.seu_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(vis.seu_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

SpatialDimPlot(vis.seu,alpha = c(0.1,0.1))

library(SeuratWrappers)
library(Banksy)

vis.seu <- RunBanksy(vis.seu,
                     lambda = 0.8, verbose = TRUE,
                     assay = "Spatial.016um", slot = "data", features = "variable",
                     k_geom = 50)

DefaultAssay(vis.seu) <- "BANKSY"
vis.seu <- RunPCA(vis.seu, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(vis.seu), npcs = 30)
vis.seu <- FindNeighbors(vis.seu, reduction = "pca.banksy", dims = 1:30)
vis.seu <- FindClusters(vis.seu, cluster.name = "banksy_cluster", resolution = 0.5)
Idents(vis.seu) <- "banksy_cluster"
p <- SpatialDimPlot(vis.seu, group.by = "banksy_cluster")
p

banksy_cells <- CellsByIdentities(vis.seu)
p <- SpatialDimPlot(vis.seu,ncol=6 ,cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
p

cortex.coordinates <- as.data.frame(read.csv("/brahms/lis/visium_hd/final_mouse/cortex-hippocampus_coordinates.csv"))
cortex <- CreateSegmentation(cortex.coordinates)

object[["cortex"]] <- Overlay(object[["slice1.008um"]], cortex)
cortex <- subset(object, cells = Cells(object[["cortex"]]))

###Integration with scRNA-seq data (deconvolution)

library(spacexr)
# sketch the cortical subset of the Visium HD dataset
cortex <- vis.seu
DefaultAssay(cortex) <- "Spatial.008um"
cortex <- FindVariableFeatures(cortex)
cortex <- SketchData(
  object = cortex,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(cortex) <- "sketch"
cortex <- ScaleData(cortex)
cortex <- RunPCA(cortex, assay = "sketch", reduction.name = "pca.cortex.sketch", verbose = T)
cortex <- FindNeighbors(cortex, reduction = "pca.cortex.sketch", dims = 1:50)
cortex <- RunUMAP(cortex, reduction = "pca.cortex.sketch", reduction.name = "umap.cortex.sketch", return.model = T, dims = 1:50, verbose = T)

# load in the reference scRNA-seq dataset
ref <- readRDS("allen_scRNAseq_ref.Rds")
Idents(ref) <- "subclass_label"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$subclass_label)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

counts_hd <- cortex[["sketch"]]$counts
cortex_cells_hd <- colnames(cortex[["sketch"]])
coords <- GetTissueCoordinates(cortex)[cortex_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 10)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
cortex <- AddMetaData(cortex, metadata = RCTD@results$results_df)
# project RCTD labels from sketched cortical cells to all cortical cells
cortex$first_type <- as.character(cortex$first_type)
cortex$first_type[is.na(cortex$first_type)] <- "Unknown"
cortex <- ProjectData(
  object = cortex,
  assay = "Spatial.008um",
  full.reduction = "pca.cortex",
  sketched.assay = "sketch",
  sketched.reduction = "pca.cortex.sketch",
  umap.model = "umap.cortex.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)
Idents(cortex) <- "full_first_type"

# now we can spatially map the location of any scRNA-seq cell type
# start with Layered (starts with L), excitatory neurons in the cortex
cells <- CellsByIdentities(cortex)
excitatory_names <- sort(grep("^L.* CTX", names(cells), value = TRUE))
excitatory_names <- names(cells)[c(1,5,16)]
p <- SpatialDimPlot(cortex, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 4)
p

plot_cell_types <- function(data, label) {
  p <- ggplot(data, aes(x = get(label), y = n, fill = full_first_type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = ifelse(n >= min_count_to_show_label, full_first_type, "")), position = position_stack(vjust = 0.5), size = 2) +
    xlab(label) +
    ylab("# of Spots") +
    ggtitle(paste0("Distribution of Cell Types across ", label)) +
    theme_minimal()
}

cell_type_banksy_counts <- cortex[[]] %>%
  dplyr::filter(full_first_type %in% excitatory_names) %>%
  dplyr::count(full_first_type, banksy_cluster)

min_count_to_show_label <- 20

p <- plot_cell_types(cell_type_banksy_counts, "banksy_cluster")
p

DefaultAssay(cortex) <- "Spatial.008um"
SpatialDimPlot(cortex,alpha = c(0.8,1),interactive = F,pt.size.factor = 1,group.by = c("full_first_type","banksy_cluster"))
ggsave(height = 7,width = 18,filename = "dimplot.png")
############
# set ID to RCTD label
Idents(cortex) <- "full_first_type"

# Visualize distribution of 4 interneuron subtypes
inhibitory_names <- c("Sst", "Pvalb", "Vip", "Lamp5")
cell_ids <- CellsByIdentities(cortex, idents = inhibitory_names)
p <- SpatialDimPlot(cortex, cells.highlight = na.omit(cell_ids), cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 4)
p
ggsave(filename = "inN.png")
# create barplot to show proportions of cell types of interest
layer_table <- table(object$full_first_type, object$layer_id)[inhibitory_names, 1:4]

neuron_props <- reshape2::melt(prop.table(layer_table), margin = 1)
ggplot(neuron_props, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cell type", y = "Proportion", fill = "Layer") +
  theme_classic()
