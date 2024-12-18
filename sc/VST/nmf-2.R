rm(list = ls())
setwd("~/rna/sc/VST")
int.seu <- qs::qread("vst.qs")
# mye.seu <- qs::qread("mye_new.seu.qs")
library(Seurat)
library(GeneNMF)
vst.seu <- subset(int.seu,celltype == "Schwannoma")
sampleObj.list <- Seurat::SplitObject(vst.seu, split.by = "orig.ident")
geneNMF.programs <- multiNMF(sampleObj.list, k=6:10)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs, nMP=8,weight.explained = 0.8,
                                        metric = "cosine",min.confidence = 0.7)
View(geneNMF.metaprograms$metaprograms.metrics)

plotMetaPrograms(geneNMF.metaprograms)
lapply(geneNMF.metaprograms$metaprograms.genes, head)
library(msigdbr)
library(fgsea)

top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(vst.seu), category = "C8")
})

library(UCell)
mp.genes <- geneNMF.metaprograms$metaprograms.genes
vst.seu <- AddModuleScore_UCell(vst.seu, features = mp.genes, ncores=4, name = "")

VlnPlot(vst.seu, features=names(mp.genes), group.by = "orig.ident",
        pt.size = 0, ncol=5)

matrix <- vst.seu@meta.data[,names(mp.genes)]
vst.seu$nmf_cluster <- apply(matrix,1,which.max)

#dimred <- scale(matrix)
dimred <- as.matrix(matrix)

colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
vst.seu@reductions[["MPsignatures"]] <- new("DimReduc",
                                        cell.embeddings = dimred,
                                        assay.used = "RNA",
                                        key = "MP_",
                                        global = FALSE)

set.seed(123)
vst.seu <- RunUMAP(vst.seu, reduction="MPsignatures", dims=1:length(vst.seu@reductions[["MPsignatures"]]),
               metric = "euclidean", reduction.name = "nmf")

library(viridis)
library(ggplot2)
FeaturePlot(vst.seu, features = names(mp.genes), reduction = "nmf", ncol=4,raster = F) &
  scale_color_viridis(option="B") &
  theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())

seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799","#D16258","#BDBF67","#924099")

seurat_colors <- MetBrewer::met.brewer("Klimt",5)

DimPlot(vst.seu,reduction = "umap",group.by = c("nmf_cluster","batch"),
        label = T,repel = T,label.size = 3,label.box = T,
        cols = seurat_colors,label.color = "grey100",
        pt.size = 0.01,raster = F) 
vst.seu = NormalizeData(vst.seu,verbose=T,assay = "RNA")
vst.seu = FindVariableFeatures(vst.seu,selection.method = "vst",nfeatures= 2000,verbose=F)
vst.seu = ScaleData(vst.seu,verbose=F)
vst.seu = RunPCA(vst.seu,verbose=F)
ElbowPlot(vst.seu,ndims = 50)
vst.seu$batch <- substr(colnames(vst.seu),1,4)
vst.seu <- harmony::RunHarmony(vst.seu, group.by.vars = "orig.ident",
                               ncores=10,plot_convergence=T,
                               reduction = "pca", assay.use = "RNA", 
                               kmeans_init_nstart=20, kmeans_init_iter_max=100,
                               reduction.save = "harmony",max_iter=30)
vst.seu = FindNeighbors(vst.seu,dims=seq(30),verbose=F,reduction = "harmony")
vst.seu <- RunUMAP(vst.seu,reduction = "harmony",dims = 1:30)
vst.seu <- FindClusters(vst.seu,resolution = 0.1)
DimPlot(vst.seu,group.by = c("seurat_clusters","batch"),raster = F,label = T)

qs::qsave(vst.seu,file="mye_new.seu.qs")
