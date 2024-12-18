###Sub cluster######
rm(list = ls())
setwd("/home/hzg/rna2/cri")
int.seu <- qs::qread("int.seu.qs")
library(Seurat)
library(dplyr)
DimPlot(int.seu)
library(harmony)
CMO.seu <- readRDS("/home/hzg/rna/scRNA_analysis/BT85CMO/int.seu.Rds")
int.seu <-merge(int.seu,CMO.seu)
int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.85&percent.mt<20& 
                    nFeature_RNA >=500&nCount_RNA>=1000)
# Standard PCA
int.seu <- int.seu %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst",nfLeatures= 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50,verbose = F)
int.seu <- RunHarmony(int.seu, group.by.vars="orig.ident", max_iter = 25)

ElbowPlot(int.seu, ndims=50, reduction="pca") 
pc.num <- 1:20

int.seu <- int.seu %>% 
  RunUMAP(reduction = "pca", dims = pc.num, verbose = F) %>% 
  FindNeighbors(reduction = "pca",k.param = 10, dims = pc.num)
DimPlot(int.seu,group.by = "seurat_clusters",label = T)+ggtitle("BT109")
FeaturePlot(int.seu,features = c("MET","eGFP"))

int.seu <- subset(int.seu,orig.ident%in%c("Crizotinib_long", "Crizotinib_short","DMSO"))


###Cluster######
setwd("/home/hzg/rna2/cri")
rm(list=ls())
gc()
pc.num <- 1:20
int.seu <- int.seu %>% 
  RunUMAP(reduction = "pca", dims = pc.num, verbose = F) %>% 
  FindNeighbors(reduction = "pca",k.param = 10, dims = pc.num)
int.seu <- int.seu %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()
DimPlot(int.seu, reduction = "umap", label = TRUE, pt.size = .1)

cell.prop <- as.data.frame(prop.table(table(Idents(int.seu), int.seu$orig.ident)))
colnames(cell.prop) <- c("cluster", "sample", "proportion")
ggplot(cell.prop, aes(sample, proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = sample(MetBrewer::met.brewer("Klimt",25),25))+
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.1, "cm"),
    legend.key.size = unit(10, "pt"),
    axis.text = element_text(size = 5, face = "bold"),
    title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_legend(title = NULL))

##Subcluster##
setwd("/home/hzg/rna2/cri")
rm(list=ls())
gc()
int.seu <-qs::qread("tumor.qs")
library(decontX)
counts <- int.seu@assays$RNA@counts
decontX_results <- decontX(counts) 
int.seu$Contamination <- decontX_results$contamination
FeaturePlot(int.seu, 
            features = 'Contamination', 
            raster=FALSE) + scale_color_viridis_c()+theme_bw()+theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
int.seu@assays$RNA@counts <- decontX_results$decontXcounts

int.seu <- int.seu %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst",nfeatures= 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50,verbose = F)
ElbowPlot(int.seu, ndims=50, reduction="pca") 
pc.num <- 1:30
int.seu <- int.seu %>% 
  RunUMAP(reduction = "pca", dims = pc.num, verbose = F) %>% 
  FindNeighbors(reduction = "pca",k.param = 10, dims = pc.num)
int.seu <- FindClusters(int.seu,resolution = 0.8) %>% identity()
DimPlot(int.seu,reduction = "umap", label = T,repel = F,label.size = 3,
        label.box = T,label.color = "grey100",pt.size = 0.1)+theme(legend.position = "right")

all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,max.cells.per.ident = 200,
                              logfc.threshold = 0.25) 
top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DimPlot(int.seu,reduction = "umap", label = T,repel = F,label.size = 3,
        label.box = T,label.color = "grey100",pt.size = 0.1)+theme(legend.position = "right")
qs::qsave(int.seu,file = "int_adj.qs")
