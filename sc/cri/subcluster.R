###Sub cluster######
rm(list = ls())
setwd("/home/hzg/rna/sc/cri")
int.seu <- qs::qread("tumor.qs")
library(Seurat)
library(dplyr)
DimPlot(int.seu)
library(harmony)
CMO.seu <- readRDS("/home/hzg/rna/scRNA_analysis/BT85CMO/int.seu.Rds")
int.seu <-merge(int.seu,CMO.seu)
int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.85&percent.mt<20& 
                    nFeature_RNA >=250&nCount_RNA>=1000)
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
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(int.seu, reduction = "umap", label = TRUE, pt.size = .1)

cell.prop <- as.data.frame(prop.table(table(Idents(int.seu), int.seu$state)))
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

cell.prop <- as.data.frame(prop.table(table(int.seu$H3K27M,int.seu$cnv_leiden)))
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
7