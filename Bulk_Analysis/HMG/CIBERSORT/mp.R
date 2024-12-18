library(Seurat)
library(harmony)
mg.seu <- Read10X("NDGBM/")
mg.seu <- CreateSeuratObject(counts=mg.seu,min.cells = 3,min.features = 200)
metadata <- read.csv("annot_Human_ND_GBM_TAM.csv")
mg.seu <- mg.seu[,metadata$cell]
mg.seu$cell <- rownames(mg.seu@meta.data)
metadata2 <- left_join(mg.seu@meta.data,metadata)
rownames(metadata2) <- metadata2$cell
mg.seu@meta.data <- metadata2
mds <- as.matrix(mg.seu@meta.data[,c("x","y")])
colnames(mds) <- c("UMAP_1","UMAP_2")
mg.seu[["umap"]] <- CreateDimReducObject(embeddings = mds, key = "UMAP_", assay = DefaultAssay(mg.seu))
Idents(mg.seu) <- mg.seu$cluster
DimPlot(mg.seu,label = T,group.by = "cluster")
mg.seu@active.ident
mg.downsample.seu <- subset(mg.seu,downsample=1000)

mg.seu <- subset(mg.seu,cluster%in%c("Mo-TAM","Monocytes","Unknown"),invert=T)
mg.seu <- NormalizeData(mg.seu)
mg.seu <- FindVariableFeatures(mg.seu, selection.method = "vst", nfeatures = 3000)
mg.seu <- ScaleData(mg.seu)
mg.seu <- RunPCA(mg.seu, verbose = FALSE)
ElbowPlot(mg.seu)
mg.seu <- RunHarmony(mg.seu, group.by.vars = "sample", max_iter = 20,
                      reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
mg.seu <- RunUMAP(mg.seu,reduction = "harmony",dims = 1:15)
mg.seu <- FindNeighbors(mg.seu,reduction = "harmony")
mg.seu <- FindClusters(mg.seu,resolution = 0.4)
DimPlot(mg.seu,label = T,group.by = "seurat_clusters")
Idents(mg.seu)<- mg.seu$cluster
seurat_colors <- c("#CC4E6A","#4A4999","#D77B5A" , "#5DA373" )
DimPlot(mg.seu,group.by = "cluster",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)
ggsave("umap_mg.png",height = 6,width = 9)
markers <- FindAllMarkers(mg.seu,only.pos = T,logfc.threshold = 0.8)
qs::qsave(mg.seu,"mg.seu")
mg.seu <- qs::qread("mg.seu")
##########
dat <- as.data.frame(mg.seu@assays$RNA$counts)
dat[24810,] <- paste0("MG",mg.seu$RNA_snn_res.0.4)
write.table(dat,file = "HMG_Count.txt",quote = F,sep = "\t")
write.table(as.data.frame(exprSet_filter),file = "HMG_mixture.txt",quote = F,sep = "\t")
