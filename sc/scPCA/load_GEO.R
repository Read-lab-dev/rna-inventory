setwd("~/rna/sc/scPCA/")
rm(list=ls())
gc()
library(Seurat)
library(dplyr)
packageVersion("Seurat")
path <- dir("./raw2/")
sample_names <- path
scRNAlist <- list()
project_name <- paste0("CU",clipr::read_clip())
# ensembl2symbol <- clusterProfiler::bitr(rownames(tmp),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
# ensembl2symbol <- ensembl2symbol[!duplicated(ensembl2symbol$ENSEMBL),]
for (i in 1:length(path)) {
  tmp <- Read10X(paste0("./raw2/",path[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = tmp,project = project_name[i])
  scRNAlist[[i]]$orig.ident <- project_name[i]
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=project_name[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)
int.seu <- subset(int.seu, subset = nFeature_RNA >=250&nCount_RNA<20000&percent.mt<25)
int.seu <- SCTransform(int.seu,vars.to.regress = "percent.mt")
int.seu = RunPCA(int.seu,verbose=F)
int.seu = harmony::RunHarmony(int.seu,group.by.vars="orig.ident",max_iter=20)
ElbowPlot(int.seu,ndims = 50)
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu = FindNeighbors(int.seu,reduction = "harmony",dims=seq(30),verbose=F)
int.seu = FindClusters(int.seu,res=1,verbose=F)
DimPlot(int.seu,reduction = "umap",label = T)

#######################
DimPlot(int.seu,reduction = "umap",label = T,group.by = "orig.ident")

FeaturePlot(int.seu,features = c("GFAP","OLIG1","OLIG2","MAP2","PTPRC","PDGFRA","EGFR","YAP1","WWTR1","TEAD1","TEAD4"),order = T,reduction = "umap")+
  DimPlot(int.seu,reduction = "umap",label = T)&NoLegend()&NoAxes()
ggsave("HGG-marker.png",height = 9,width = 16)

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
write.csv(cellmarker,file = "cellmarkerHGG.csv")
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
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",10),10)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)
ggsave("UMAP_HGG.png",dpi = 200)

DimPlot(int.seu, group.by = "treat",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = sample(MetBrewer::met.brewer("Klimt",4),2),label.color = "grey100",
        pt.size = 0.1)

ggsave("UMAP_sham.png",dpi = 200)

DimPlot(int.seu,label = T)

qs::qsave(int.seu,file = "./scPCA_LGG.qs")

########Module#########
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gghighlight)
library(MetBrewer)
library(RColorBrewer)
library(tibble)
lineage <- read.csv("../../BT109VP/NMF/module4.csv",header = T,na.strings = "")
gsc.seu <- subset(int.seu,celltype=="Neoplastic")
gsc.seu <- subset(gsc.seu,cnv.score>0.02)

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
lineage.data$YAP1 <- GetAssayData(gsc.seu,assay = "SCT")["TEAD1",]
seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=orig.ident))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = MetBrewer:::met.brewer("Klimt",19))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))
ggsave("lineage_cnv.png",dpi = 200)
ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=lineage.class))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors[c(2,4,5)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))
ggsave("lineage_class_cnv.png",dpi = 200)

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=YAP1))+
  geom_point(alpha=1,size=0.8)+scale_color_gradient2(low = "black",mid = "grey",high = "Purple",midpoint = 0)+
  theme_few()+xlab("AC-like <--------------> NPC-like")
