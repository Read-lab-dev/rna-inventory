.libPaths(c("/home/hzg/R/x86_64-pc-linux-gnu-library/SeuratV5/","/home/hzg/R/x86_64-pc-linux-gnu-library/4.2"))
detach("package:Seurat", unload = T)
options(future.globals.maxSize= 891289600)
library(Seurat)
setwd("/home/hzg/rna/visium/spatial/xenium/FFPE/")
xenium.obj <- LoadXenium("/home/hzg/rna/visium/spatial/xenium/FFPE")

xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

ImageDimPlot(xenium.obj, fov = "fov", molecules  = c("CD4", "APOE", "AQP4", "OLIG1","CD8A"), nmols = 20000)

head(GetTissueCoordinates(xenium.obj[["fov"]][["centroids"]]))
cropped.coords <- Crop(xenium.obj[["fov"]], x = c(3000, 6800), y = c(0, 3000), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("CD4", "APOE", "AQP4", "OLIG1","CD8A"),alpha = 0.5)

xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:8)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:8)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.2)

DimPlot(xenium.obj)


xenium.obj <- qs::qread("xenium.obj.qs")
FeaturePlot(xenium.obj,c("APOE", "AQP4", "OLIG1","EGFR","CD4","CD8A","PDGFRA","CD14","PTPRC","TOP2A"))
ImageDimPlot(xenium.obj, fov = "zoom",cols = "polychrome", size = 0.6)
allmarkers <- FindAllMarkers(xenium.obj,logfc.threshold = 0.5)
library(dplyr)
top10.markers <- allmarkers%>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
ImageFeaturePlot(xenium.obj,fov = "zoom",c("APOE", "AQP4", "OLIG1","EGFR","CD4","CD8A","PDGFRA","CD14","PTPRC","TOP2A"),blend=F)

library(spacexr)
query.counts <- GetAssayData(xenium.obj, assay = "Xenium", slot = "counts")[, Cells(xenium.obj[["zoom"]])]
coords <- GetTissueCoordinates(xenium.obj[["zoom"]], which = "centroids")
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
aa <- intersect(toupper(gene_ids$SYMBOL),rownames(xenium.obj))
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
xenium.obj$predicted.celltype <- annotations
keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
xenium.obj <- subset(xenium.obj, cells = keep.cells)

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "zoom", group.by = "predicted.celltype",niches.k = 5, neighbors.k = 30)

qs::qsave(xenium.obj,file = "xenium.obj.qs")
qs::qsave(myRCTD,file = "xenium_RTCD.obj.qs")


head(GetTissueCoordinates(xenium.obj[["fov"]][["centroids"]]))
cropped.coords <- Crop(xenium.obj[["fov"]], x = c(3000, 6800), y = c(0, 3000), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
ImageFeaturePlot(xenium.obj,fov = "zoom",c("APOE", "AQP4", "OLIG1","EGFR","CD4","CD8A","PDGFRA","CD14","PTPRC","TOP2A"))
ImageFeaturePlot(xenium.obj,fov = "zoom",features =c("DCN"))
