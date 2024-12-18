.libPaths(c("/home/hzg/R/x86_64-pc-linux-gnu-library/SeuratV5/","/home/hzg/R/x86_64-pc-linux-gnu-library/4.4"))
detach("package:Seurat", unload = T)
library(Seurat)
use_condaenv("r-reticulate",required = T)
library(reticulate)
# library(SeuratDisk)
# Convert("adata_16samples_DMG.h5ad", dest = "h5seurat", overwrite = TRUE)
liu.obj <- LoadH5Seurat("adata_16samples_DMG.h5seurat")
liu.obj <- RunPCA(liu.obj, features = rownames(liu.obj))

liu.obj <- FindNeighbors(liu.obj, dims = 1:10)
liu.obj <- FindClusters(liu.obj, resolution = 0.5)

liu.obj <- RunUMAP(liu.obj, dims = 1:10)
DimPlot(liu.obj, reduction = "umap")

spatial_coords <- data.frame(
  x = liu.obj$X,
  y = liu.obj$Y,
  row.names = colnames(liu.obj)  # Ensure row names match the cell barcodes
)

spatial_image <- CreateSeuratImage(
  spatial_coords, 
  dims = c("x", "y"),
  image.name = "example"
)

liu.obj <- addimage(liu.obj, metadata = spatial_coords)

ImageDimPlot(liu.obj,group.by = "name")

ImageDimPlot(liu.obj, fov = "fov", molecules  = c("CD4", "APOE", "AQP4", "OLIG1","CD8A"), nmols = 20000)

head(GetTissueCoordinates(liu.obj[["fov"]][["centroids"]]))
cropped.coords <- Crop(liu.obj[["fov"]], x = c(3000, 6800), y = c(0, 3000), coords = "plot")
liu.obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(liu.obj[["zoom"]]) <- "segmentation"
ImageDimPlot(liu.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("CD4", "APOE", "AQP4", "OLIG1","CD8A"),alpha = 0.5)

liu.obj <- SCTransform(liu.obj, assay = "Xenium")
liu.obj <- RunPCA(liu.obj, npcs = 30, features = rownames(liu.obj))
liu.obj <- RunUMAP(liu.obj, dims = 1:8)
liu.obj <- FindNeighbors(liu.obj, reduction = "pca", dims = 1:8)
liu.obj <- FindClusters(liu.obj, resolution = 0.2)

DimPlot(liu.obj)


liu.obj <- qs::qread("liu.obj.qs")
FeaturePlot(liu.obj,c("APOE", "AQP4", "OLIG1","EGFR","CD4","CD8A","PDGFRA","CD14","PTPRC","TOP2A"))
ImageDimPlot(liu.obj, fov = "zoom",cols = "polychrome", size = 0.6)
allmarkers <- FindAllMarkers(liu.obj,logfc.threshold = 0.5)
library(dplyr)
top10.markers <- allmarkers%>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
ImageFeaturePlot(liu.obj,fov = "zoom",c("APOE", "AQP4", "OLIG1","EGFR","CD4","CD8A","PDGFRA","CD14","PTPRC","TOP2A"),blend=T)

library(spacexr)
query.counts <- GetAssayData(liu.obj, assay = "Xenium", slot = "counts")[, Cells(liu.obj[["zoom"]])]
coords <- GetTissueCoordinates(liu.obj[["zoom"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

allen.cortex.ref <- readRDS("./allen_cortex.rds")
allen.cortex.ref <- UpdateSeuratObject(allen.cortex.ref)

library(org.Mm.eg.db)
# allen.corted.ref can be downloaded here:
# https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1
allen.cortex.ref <- readRDS("/brahms/shared/vignette-data/allen_cortex.rds")
allen.cortex.ref <- UpdateSeuratObject(allen.cortex.ref)
gene_ids<-AnnotationDbi::select(org.Mm.eg.db, keys=as.character(rownames(allen.cortex.ref)), 
                                columns="SYMBOL", #目标格式
                                keytype="ALIAS") #目前的格式
gene_ids
aa <- intersect(toupper(gene_ids$SYMBOL),rownames(liu.obj))
aa <- stringr::str_to_title(aa)
allen.cortex.ref <- subset(allen.cortex.ref,features=aa)
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
allen.cortex.ref <- RenameGenesSeurat(allen.cortex.ref,newnames = toupper(rownames(allen.cortex.ref)))
Idents(allen.cortex.ref) <- "subclass"
# remove CR cells because there aren't enough of them for annotation
allen.cortex.ref <- subset(allen.cortex.ref, subset = subclass != "CR")
counts <- GetAssayData(allen.cortex.ref, assay = "RNA", slot = "counts")
cluster <- as.factor(allen.cortex.ref$subclass)
names(cluster) <- colnames(allen.cortex.ref)
nUMI <- allen.cortex.ref$nCount_RNA
names(nUMI) <- colnames(allen.cortex.ref)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)
# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 1)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
liu.obj$predicted.celltype <- annotations
keep.cells <- Cells(liu.obj)[!is.na(liu.obj$predicted.celltype)]
liu.obj <- subset(liu.obj, cells = keep.cells)

liu.obj <- BuildNicheAssay(object = liu.obj, fov = "zoom", group.by = "predicted.celltype",niches.k = 5, neighbors.k = 30)

qs::qsave(liu.obj,file = "liu.obj.qs")
qs::qsave(myRCTD,file = "xenium_RTCD.obj.qs")


head(GetTissueCoordinates(liu.obj[["fov"]][["centroids"]]))
cropped.coords <- Crop(liu.obj[["fov"]], x = c(3000, 6800), y = c(0, 3000), coords = "plot")
liu.obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(liu.obj[["zoom"]]) <- "segmentation"
ImageFeaturePlot(liu.obj,fov = "zoom",c("APOE", "AQP4", "OLIG1","EGFR","CD4","CD8A","PDGFRA","CD14","PTPRC","TOP2A"))
ImageFeaturePlot(liu.obj,fov = "zoom",features =c("DCN"))
