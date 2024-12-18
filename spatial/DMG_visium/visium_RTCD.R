rm(list=ls())
gc()
setwd("~/rna/spatial/DMG_visium")
library(Seurat)
vis.seu <- qs::qread("./seurat/DMG.merge.qs")
library(spacexr)
library(dplyr)
sc <- qs::qread("~/rna/BT109VP/gsc.seu.qs")
sc <- subset(sc,orig.ident=="Eng_DMSO")
### Load in/preprocess your data, this might vary based on your file type
#refdir <- system.file("extdata",'Reference/Vignette',package = 'RCTD') # directory for the reference
counts <- as.matrix(sc@assays$RNA@counts)#read.csv(file.path(refdir,"dge.csv")) # load in counts matrix
#rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
meta_data <- sc@meta.data #read.csv(file.path(refdir,"meta_data.csv")) # load in meta_data (barcodes, clusters, and nUMI)
cell_types <- meta_data$celltype; names(cell_types) <- rownames(meta_data)#meta_data$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- apply(counts, 2, sum); names(nUMI) <- rownames(meta_data)#meta_data$barcode # create nUMI named list
mode(counts) <- "integer";mode(nUMI) <- "integer"
### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)

spatialcounts <- as.matrix(vis.seu@assays$Spatial$counts)
coords <- GetTissueCoordinates(vis.seu,scale=NULL)
coords <- bind_rows(lapply(Images(vis.seu), function(x) GetTissueCoordinates(vis.seu, image = x, scale = NULL)))

UMI <- vis.seu@meta.data$nCount_Spatial
names(UMI) <- rownames(vis.seu@meta.data)
puck <- SpatialRNA(coords = coords, counts=spatialcounts, nUMI=UMI)

# Sys.setenv("OPENBLAS_NUM_THREADS"=8)
myRCTD <- create.RCTD(puck, reference, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet') #doublet_mode = 'doublet'
vis.seu <- AddMetaData(vis.seu, metadata = myRCTD@results$results_df)
cell_to_use <- colnames(vis.seu)[vis.seu$first_type!="NA"]
vis.seu <- subset(vis.seu,cells =cell_to_use)
seurat_colors <- as.character(MetBrewer::met.brewer("VanGogh3", 7))
seurat_colors <- c("yellow", "skyblue", "red", "green")
SpatialDimPlot(vis.seu,group.by = c("first_type","second_type"),ncol = 5)&scale_fill_manual(values = seurat_colors)
ggsave(filename = "RCTD-DMG.png",height = 9,width = 16)
qs::qsave(vis.seu,file = "DMG.RTCD.qs")
