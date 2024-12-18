rm(list = ls())
setwd("/home/hzg/backup/sc_analysis/BT85CMO/")
library(Seurat)
library(ggplot2)
library(dplyr)
tab <- read.csv("/home/hzg/rna/seq/BT109sh/cmo/outs/multi/multiplexing_analysis/assignment_confidence_table.csv")
seurat_data <- Read10X(data.dir = "/home/hzg/rna/seq/BT109sh/cmo/outs/multi/count/raw_feature_bc_matrix/")
int.seu <- CreateSeuratObject(counts = seurat_data$`Gene Expression`,
                                   min.features = 250,
                                   min.cells = 3)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)

VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)

FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")+
  FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.85&percent.mt<20& 
                    nFeature_RNA >=250&nCount_RNA>=1000)
# Standard PCA
int.seu <- int.seu %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst",nfLeatures= 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50,verbose = F)

ElbowPlot(int.seu, ndims=50, reduction="pca") 
pc.num <- 1:30

int.seu <- int.seu %>% 
  RunUMAP(reduction = "pca", dims = pc.num, verbose = F)

int.seu <- int.seu %>% FindNeighbors(reduction = "pca") %>% 
  FindClusters(resolution = 0.3)
DimPlot(int.seu, reduction = "umap", label = TRUE, pt.size = .1)
FeaturePlot(int.seu,features = "ACTA2")

all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 1) 

top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  
  cluster <- paste(tmp$gene,collapse = ",")
  
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker.csv")

int.seu <- subset(int.seu,idents=c(0:4,7:13))

new.cluster.ids <- readLines("clipboard")##Copy from your handmade-Celltype
Idents(int.seu) <- int.seu$seurat_clusters
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
DimPlot(int.seu, reduction = "umap", label = F,repel = T,label.size = 4,
        label.box = T,label.color = "black",raster = T,raster.dpi = c(1000,1000),
        pt.size = 2.5,cols = rainbow(11,rev = T))+theme(legend.position = "right")

cell.prop<-as.data.frame(prop.table(table(Idents(int.seu), int.seu$orig.ident)))

colnames(cell.prop)<- c("cluster","sample","proportion")

ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  scale_fill_manual(values = rainbow(11,rev = T))+
  geom_bar(stat="identity",position="fill")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))

saveRDS(int.seu,file = "int.seu.Rds")

tab <- left_join(int.seu@meta.data,tab)
tab$Assignment[is.na(tab$Assignment)] <- "Unassigned"
int.seu$assignment <- tab$Assignment

#############OWN###########
# features <- list(
#   Radialglia = c("SLC1A3","NES"),
#   NPC = c("EOMES","MKI67"),
#   OPC    = c("PDGFRA","CSPG4"),
#   Oligo      = c("OLIG1","OLIG2"),
#   Astrocyte  = c("SLC1A2","S100B","ALDH1L1"),
#   MatureNeuron       = c("RBFOX3","MAP2","DLG4"),
#   ImatureNeuron       = c("TBR1","DCX","NEUROD1"),
#   GABANeuron = c("SLC6A1","GAD2","GAD1"),
#   GluNeuron=c("SLC17A7","SLC17A6","GRIN1","GRIN2B"),
#   DopaNeuron= c("TH","SLC6A3","FOXA2","KCNJ6","NR4A2","LZMX1B"),
#   SerotonergicNeuron=c("TPH1","SLC6A4","FEV"),
#   CholinergicNeuron=c("CHAT","SLC18A3","ACHE")
# )
# #Sloan Lab
# features <- list(
#   Astrocyte = c("HES1","GFAP","AGT","GLUL"),
#   Neuroepithelial = c("KLF4","VIM","VCAN"),
#   Neuron    = c("RBFOX3","DCX","MAP2"),
#   pNPC      = c("PAX6","SOX2","GLI3","TOP2A"),
#   EarlyNeuroectodermProgenitors  = c("PCP4","NEFL","NEFM"),
#   Fib       = c("COL1A1","PDGFRA","DCN"),
#   GliomaStemCell= c("EGFR","CD44","PROM1","FUT4","NEU1","NES","eGFP"))
DotPlot(int.seu,features = features)

