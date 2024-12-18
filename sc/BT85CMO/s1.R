rm(list = ls())
setwd("/home/hzg/rna/sc/BT85CMO/")
library(Seurat)
library(ggplot2)
library(dplyr)
samples <- list.files("rawdata")
seurat_list <- list()

for (sample in samples) {
  data.path <- paste0("./rawdata/", sample,"/count/sample_filtered_feature_bc_matrix")
  seurat_data <- Read10X(data.dir = data.path)
  seurat_obj <- CreateSeuratObject(counts = seurat_data$`Gene Expression`,
                                   project = sample,
                                   min.features = 200,
                                   min.cells = 3)
  seurat_list <- append(seurat_list, seurat_obj)
}
int.seu <- merge(seurat_list[[1]],y = seurat_list[-1],add.cell.ids = samples)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)

VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)

FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")+
  FeatureScatter(int.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.85&percent.mt<40)
int.seu <-JoinLayers(int.seu)
int.seu <- subset(int.seu,orig.ident%in%c("ORG_DMSO24","ORG_VP14","ORG_VP24"))
# Standard PCA
int.seu <- int.seu %>% 
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst",nfLeatures= 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50,verbose = F)

ElbowPlot(int.seu, ndims=50, reduction="pca") 
pc.num <- 1:20

int.seu <- int.seu %>% 
  RunUMAP(reduction = "pca", dims = pc.num, verbose = F) %>% 
  FindNeighbors(reduction = "pca", k.param = 10, dims = pc.num)

library(clustree)
obj <- FindClusters(int.seu, resolution = seq(0.3,1.0,by=0.1))
clustree(obj)
rm(obj)

int.seu <- int.seu %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()
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
write.csv(cellmarker,file = "cellmarker_new.csv")

int.seu <- subset(int.seu,idents=c(0:4,7:13))

new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
Idents(int.seu) <- int.seu$seurat_clusters
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
seurat_colour <- sample(MetBrewer::met.brewer("Klimt",10),5)
p1= DimPlot(int.seu, reduction = "umap", label = F,repel = T,label.size = 4,
        label.box = T,label.color = "black",raster = F,raster.dpi = c(1000,1000),
        pt.size = 0.5,cols = seurat_colour)+theme(legend.position = "right")

int.seu$orig.ident <- gsub("ORG","Eng",int.seu$orig.ident)
int.seu$orig.ident <- gsub("24","24h",int.seu$orig.ident)
int.seu$orig.ident <- gsub("14","14d",int.seu$orig.ident)
int.seu$orig.ident <- factor(int.seu$orig.ident,levels = c("Eng_DMSO24h","Eng_VP24h","Eng_VP14d"))

cell.prop<-as.data.frame(prop.table(table(Idents(int.seu),int.seu$orig.ident)))

colnames(cell.prop)<- c("cluster","sample","proportion")

p2 <-ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  scale_fill_manual(values = seurat_colour)+
  geom_bar(stat="identity",position="fill")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))

p1+p2+patchwork::plot_layout(widths = c(2,1))
ggsave("porp.pdf",height = 4,width = 8)
saveRDS(int.seu,file = "int.seu.Rds")
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

