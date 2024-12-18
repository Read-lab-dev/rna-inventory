# install.packages("devtools")
libpath <- .libPaths()
.libPaths(c("/home/hzg/R/x86_64-pc-linux-gnu-library/SeuratV5/",libpath))
rm(list = ls())
library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)
setwd("~/rna/visium")

sc <- qs::qread("~/rna/BT109VP/gsc.seu.qs")
### Load in/preprocess your data, this might vary based on your file type
#refdir <- system.file("extdata",'Reference/Vignette',package = 'RCTD') # directory for the reference
counts <- as.matrix(sc@assays$RNA@counts)#read.csv(file.path(refdir,"dge.csv")) # load in counts matrix
#rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
meta_data <- sc@meta.data #read.csv(file.path(refdir,"meta_data.csv")) # load in meta_data (barcodes, clusters, and nUMI)
cell_types <- meta_data$seurat_clusters; names(cell_types) <- rownames(meta_data)#meta_data$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- apply(counts, 2, sum); names(nUMI) <- rownames(meta_data)#meta_data$barcode # create nUMI named list
mode(counts) <- "integer";mode(nUMI) <- "integer"
### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)#> Warning in Reference(counts, cell_types, nUMI): Reference: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this
#>             is intended, there is no problem.
## Examine reference object (optional)
print(dim(reference@counts)) #observe Digital Gene Expression matrix
#> [1] 384 475
table(reference@cell_types) 
saveRDS(reference,'ref.rds')

dmg.seu <- qs::qread("./Org/seurat/DMG1.qs")
spatialcounts <- as.matrix(dmg.seu@assays$Spatial$counts)
coords <- GetTissueCoordinates(dmg.seu,scale=NULL)
UMI <- dmg.seu@meta.data$nCount_Spatial
names(UMI) <- rownames(dmg.seu@meta.data)
####Xenium##
dmg.seu <- qs::qread("/home/hzg/rna/spatial/saskia_GBM_xenium/xenium.obj.qs")
spatialcounts <- as.matrix(dmg.seu@assays$Xenium$counts)
coords <- GetTissueCoordinates(dmg.seu,scale=NULL)
coords <- tibble::column_to_rownames(coords,var = "cell")
UMI <- dmg.seu@meta.data$nCount_Xenium
names(UMI) <- rownames(dmg.seu@meta.data)
########
puck <- SpatialRNA(coords = coords, counts=spatialcounts, nUMI=UMI)

# Sys.setenv("OPENBLAS_NUM_THREADS"=8)
myRCTD <- create.RCTD(puck, reference, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet') #doublet_mode = 'doublet'
dmg.seu <- AddMetaData(dmg.seu, metadata = myRCTD@results$results_df)
seurat_colors <- as.character(MetBrewer::met.brewer("Klimt", 8))

ImageDimPlot(dmg.seu,group.by = "first_type")
ImageDimPlot(dmg.seu, fov = "fov", molecules  = c("CD4", "APOE", "AQP4", "OLIG1","CD8A"), nmols = 20000)
ImageFeaturePlot(xenium.obj, fov = "fov",features = c("DCN","CD14"))



results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names

 ## you may change this to a more accessible directory on your computer.
for(i in 1:length(cell_type_names)){
  plot_puck_continuous(myRCTD@spatialRNA, colnames(dmg.seu), norm_weights[,cell_type_names[i]], title =cell_type_names[i], size=0.01,small_point=F)
  ggsave(paste(cell_type_names[i],'_weights1.jpg', sep=''), height=10, width=10, units='in', dpi=300)
}

# make the plots 
# Plots the confident weights for each cell type as in full_mode (saved as 
# 'results/cell_type_weights_unthreshold.pdf')
plot_weights(cell_type_names, puck, "./", norm_weights,resultsdir=) 
# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
# Plots the weights for each cell type as in doublet_mode. (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                     results$results_df) 

# Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

qs::qsave(myRCTD,file = "RCTD_GBM.qs")

xenium.obj <- qs::qread("./xenium.obj.qs")

Idents(dmg.seu) <- dmg.seu$first_type
npc <-WhichCells(dmg.seu,idents="Neu-like")
p1 <- ImageDimPlot(dmg.seu, group.by = "second_type")

p2 <- ImageDimPlot(dmg.seu, group.by = "first_type")
p1+p2

markers <- FindMarkers(dmg.seu, ident.1 = "Neu-like")

ImageFeaturePlot(dmg.seu,features = "DCN")

save(dmg.seu@meta.data,file="RTCD_METADATA.Rdata")
