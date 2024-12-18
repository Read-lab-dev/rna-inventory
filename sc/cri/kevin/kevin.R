setwd("~/rna2/cri/kevin")

kevin.seu <- LoadH5Seurat("analysis_scRNAseq_tumor_counts.h5seurat")
kevin.seu@meta.data <- read.csv("tumor_metadata.csv",row.names = 1)
rownames(kevin.seu@meta.data) <- colnames(kevin.seu)
kevin.seu <- subset(kevin.seu,subset = case_barcode%in%c("SM006","SM011","SM012","SM017","SM018"))
tumor.seu <- qs::qread("../tumor.qs")
int.seu <- merge(tumor.seu,kevin.seu)
int.seu$cellstate <- c(tumor.seu$state,kevin.seu$cell_state)
int.seu$batch <- c(tumor.seu$orig.ident,kevin.seu$case_barcode)

library(harmony)
int.seu <- int.seu %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst",nfLeatures= 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50,verbose = F)
ElbowPlot(int.seu, ndims=50, reduction="pca") 
pc.num <- 1:20

int.seu <- RunHarmony(int.seu, group.by.vars="batch", max_iter = 25)

int.seu <- int.seu %>% 
  RunUMAP(reduction = "harmony", dims = pc.num, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",k.param = 10, dims = pc.num)
int.seu <- FindClusters(int.seu,resolution = 0.8) %>% identity()
DimPlot(int.seu,reduction = "umap", label = T,repel = F,label.size = 3,group.by = "cellstate",raster = F,
        label.box = T,label.color = "grey100",pt.size = 0.01)+theme(legend.position = "right")
