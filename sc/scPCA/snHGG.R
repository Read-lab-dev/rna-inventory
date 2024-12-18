##Subcluster##
.libPaths()
setwd("~/rna/sc/scPCA/")
rm(list=ls())
gc()
library(Seurat)
library(dplyr)
library(ggplot2)
library(gghighlight)
library(ggthemes)
packageVersion("Seurat")
path <- grep("SCPC",dir("./HGG-nuclus/"),value=T)
sample_names <- path
scRNAlist <- list()

ensembl2symbol <- clusterProfiler::bitr(rownames(tmp),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
ensembl2symbol <- ensembl2symbol[!duplicated(ensembl2symbol$ENSEMBL),]
ensembl2symbol <- ensembl2symbol[!duplicated(ensembl2symbol$SYMBOL),]

for (i in 1:length(path)) {
  tmp <- readRDS(paste0("./HGG-nuclus/",path[i],"/",grep("processed",dir(paste0("./HGG-nuclus/",path[i])),value = T)))
  ct <- assay(tmp)[ensembl2symbol$ENSEMBL,]
  rownames(ct) <- ensembl2symbol$SYMBOL
  scRNAlist[[i]] <- CreateSeuratObject(counts = ct,meta.data =as.data.frame(colData(tmp)))
  scRNAlist[[i]]$orig.ident <- sample_names[i]
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)

int.seu <- NormalizeData(int.seu)
int.seu <- FindVariableFeatures(int.seu,nfeatures=3000)
int.seu <- ScaleData(int.seu)
int.seu <- RunPCA(int.seu)
ElbowPlot(int.seu,ndims = 50)
# int.seu <- RunUMAP(int.seu,reduction = "pca",dims = 1:30)
# int.seu = FindNeighbors(int.seu,reduction = "pca",dims=seq(30),verbose=F)
# int.seu = FindClusters(int.seu,res=0.5,verbose=F)

int.seu = harmony::RunHarmony(int.seu,group.by.vars="orig.ident",max_iter=20)
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu = FindNeighbors(int.seu,reduction = "harmony",dims=seq(30),verbose=F)
int.seu = FindClusters(int.seu,res=0.3,verbose=F)
DimPlot(int.seu,reduction = "umap",label = T)

#######################
DimPlot(int.seu,reduction = "umap",label = T,group.by = "orig.ident")

FeaturePlot(int.seu,features = c("GFAP","OLIG1","OLIG2","MAP2","PTPRC","PDGFRA","EGFR","YAP1","WWTR1","TEAD1","TEAD4"),order = T,reduction = "umap")+
  DimPlot(int.seu,reduction = "umap",label = T)

FeaturePlot(int.seu,features = c("PTPRC","EPCAM","PECAM1","MME","CD3G","CD3E","CD79A"),order = T)

all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 1) 
top10.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  
  cluster <- paste(tmp$gene,collapse = ",")
  
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker13.csv")
features <- list(
  OPC    = c("PDGFRA","CSPG4"),
  Oligo      = c("OLIG1","OLIG2","MBP"),
  Astrocyte  = c("S100B","ALDH1L1","AQP4","GFAP"),
  Neruoepithial= c("SOX2","NOTCH1","HES1"),
  Radialglia = c("SLC1A3","NES"),
  MatureN = c("RBFOX3","MAP2","DLG4"),
  Excitatory=c("SLC17A7","RORB","PCP4"),
  Inhibitory=c("GAD1","GAD2","LHX6","PVALB","LAMP5","CXCL14","CALB2"),
  MicroGlia =c("PTPRC","ITGAM","CD14"),
  Pericyte= c("PDGFRB","PECAM1"),
  Endo = c("CD34"))

DotPlot(int.seu,features = features)

feature_neu <- list(
  Excitatory=c("SLC17A7","RORB","PCP4"),
  Inhibitory=c("GAD1","GAD2","LHX6","PVALB","LAMP5","CXCL14","CALB2")
)
DotPlot(int.seu,features = feature_neu,split.by = "celltype")

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",20),7)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)
ggsave("UMAP_snHGG.png",dpi = 200)

DimPlot(int.seu, group.by = "treat",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = sample(MetBrewer::met.brewer("Klimt",4),2),label.color = "grey100",
        pt.size = 0.1)

ggsave("UMAP_sham.png",dpi = 200)

DimPlot(int.seu,label = T,split.by = "orig.ident")

qs::qsave(int.seu,file = "./snHGG.qs")



########Module#########
lineage <- read.csv("../../BT109VP/NMF/module4.csv",header = T,na.strings = "")
gsc.seu <- subset(int.seu,orig.ident%in%c("SCPCS000635","SCPCS000636","SCPCS000637","SCPCS000642"))
gsc.seu <- subset(int.seu,celltype=="Neoplastic")

gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$NPC)),nbin = 30,ctrl = 100,name = "NPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$AC)),nbin = 30,ctrl = 100,name = "AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$OPC)),nbin = 30,ctrl = 100,name = "OPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$cyc)),nbin = 30,ctrl = 100,name = "cycscore")
lineage.score <- pmax(gsc.seu$NPC1,gsc.seu$AC1)
lineage.score.plot <- lineage.score
lineage.class <- ifelse(gsc.seu$NPC1>gsc.seu$AC1,"NPC","AC")
lineage.data <- data.frame(lineage.score,lineage.score.plot,lineage.class)
lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"] <- -lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"]
lineage.data$lineage.class <- ifelse(gsc.seu$OPC1>pmax(gsc.seu$AC1,gsc.seu$NPC1),"OPC",lineage.class)
lineage.data$stemness <- gsc.seu$OPC1-lineage.data$lineage.score
lineage.data <- data.frame(lineage.data,gsc.seu@meta.data)
seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=orig.ident))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = MetBrewer:::met.brewer("Klimt",15))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))
ggsave("2dim-snHGG.png",dpi = 200)

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=lineage.class))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors[c(2,4,6)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))
ggsave("2dim-snHGG-cellstate.png",dpi = 200)
