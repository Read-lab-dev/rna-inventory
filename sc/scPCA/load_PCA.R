##Subcluster##
.libPaths()
setwd("~/rna/sc/scPCA/")
rm(list=ls())
gc()
library(Seurat)
library(dplyr)
packageVersion("Seurat")
path <- grep("SCPC",dir("./raw/"),value=T)
sample_names <- path
scRNAlist <- list()

ensembl2symbol <- clusterProfiler::bitr(rownames(tmp),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
ensembl2symbol <- ensembl2symbol[!duplicated(ensembl2symbol$ENSEMBL),]

for (i in 1:length(path)) {
  tmp <- readRDS(paste0("./raw/",path[i],"/",gsub("CS","CL",path[i]),"_processed.rds"))
  ct <- assay(tmp)[ensembl2symbol$ENSEMBL,]
  rownames(ct) <- ensembl2symbol$SYMBOL
  scRNAlist[[i]] <- CreateSeuratObject(counts = ct,meta.data =as.data.frame(colData(tmp)))
  scRNAlist[[i]]$orig.ident <- path[i]
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)

int.seu <- SCTransform(int.seu,vars.to.regress = c("subsets_mito_percent"))
int.seu = RunPCA(int.seu,verbose=F)
ElbowPlot(int.seu,ndims = 50)
int.seu <- RunUMAP(int.seu,reduction = "pca",dims = 1:30)
int.seu = FindNeighbors(int.seu,dims=seq(30),verbose=F)
int.seu = FindClusters(int.seu,res=0.5,verbose=F)
DimPlot(int.seu,reduction = "umap",label = T)

#######################
DimPlot(int.seu,reduction = "umap",label = T,group.by = "orig.ident")

FeaturePlot(int.seu,features = c("GFAP","OLIG1","OLIG2","MAP2","PTPRC","PDGFRA","EGFR","YAP1","WWTR1","TEAD1","TEAD4"),order = T,reduction = "umap")+
  DimPlot(int.seu,reduction = "umap",label = T,group.by = "orig.ident")

all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 2) 
top10.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  
  cluster <- paste(tmp$gene,collapse = ",")
  
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker22.csv")
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
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",20),20)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)
ggsave("UMAP.png",dpi = 200)

DimPlot(int.seu, group.by = "treat",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = sample(MetBrewer::met.brewer("Klimt",4),2),label.color = "grey100",
        pt.size = 0.1)

ggsave("UMAP_sham.png",dpi = 200)

DimPlot(int.seu,label = T)

qs::qsave(int.seu,file = "./scPCA.qs")

int.seu <- subset(int.seu,orig.ident%in%c("SCPCS000001","SCPCS000004","SCPCS000005","SCPCS000011","SCPCS000016","SCPCS000021"))


########Module#########
lineage <- read.csv("../../BT109VP/NMF/module4.csv",header = T,na.strings = "")
gsc.seu <- subset(int.seu,seurat_clusters%in%c(2,7,6))

gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$NPC)),nbin = 25,ctrl = 100,name = "NPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$AC)),nbin = 25,ctrl = 100,name = "AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$OPC)),nbin = 25,ctrl = 100,name = "OPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$cyc)),nbin = 25,ctrl = 100,name = "cycscore")
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
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = MetBrewer:::met.brewer("Klimt",5))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))
