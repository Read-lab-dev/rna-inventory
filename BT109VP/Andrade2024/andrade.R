rm(list = ls())
gc()
options('future.globals.maxSize'=1024^4)
library(Seurat)
library(dplyr)
library(ggplot2)
setwd("~/rna/BT109VP/Andrade2024")
path <- dir("./raw/")
sample_names <- path
scRNAlist <- NULL
for (i in 1:length(path)) {
  tmp <- Read10X(paste0("./raw/",path[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = tmp,
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^mt-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)
int.seu <- subset(int.seu, subset = nFeature_RNA >=250&percent.mt<20)
int.seu = NormalizeData(int.seu,verbose=F,assay = "RNA")
int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=F)
int.seu = ScaleData(int.seu,verbose=F)
int.seu = RunPCA(int.seu,verbose=F)
ElbowPlot(int.seu,ndims = 50)
int.seu <- harmony::RunHarmony(int.seu, group.by.vars = "orig.ident", 
                               reduction = "pca", assay.use = "RNA", 
                               reduction.save = "harmony",max_iter=20)
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu = FindNeighbors(int.seu,dims=seq(30),verbose=F,reduction = "harmony")
int.seu <- FindClusters(int.seu,resolution = 0.1)
int.seu <- JoinLayers(int.seu)
qs::qsave(int.seu,"sch.seu.qs")
DimPlot(int.seu,reduction = "umap",label = T)
DimPlot(int.seu,reduction = "umap",label = T,group.by = c("treat","technique"))

int.seu$treat <- ifelse(int.seu$orig.ident%in%c("GSM7747175","GSM7747176"),"denovo","engrafted")
int.seu$technique <- ifelse(int.seu$orig.ident%in%c("GSM7747177","GSM7747178"),"fresh","frozen")

all.markers <- FindAllMarkers(int.seu, only.pos = T,logfc.threshold = 2, min.pct = 0.25)
top10.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarkers.csv")

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype

num <- length(levels(Idents(int.seu)))
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",num),num)
features <- list(
  mo  = c("NF2","ADGRG6","SOX10","FOXD3"),
  m1  = c("YAP1","WWTR1","AXL","MERTK","TYRO3","MET"),
  m3  = c("SOX2","NES","CD3G","MS4A1","AIF1","TMEM119"))
DimPlot(int.seu,reduction = "umap",
        label = T,repel = T,label.size = 4,label.box = T,
        cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)+
  DotPlot(int.seu,features = features)
DotPlot(int.seu,features = stringr::str_to_title(unlist(features)))
library(ggplot2)
FeaturePlot(int.seu,features = stringr::str_to_title(unlist(features)),order = T,
            reduction = "umap",raster = F)&scale_colour_viridis_c(option = "B")
ggsave("mouse_featureplot.png",height = 12,width = 16)

features <-c("EPCAM","PECAM1","MME","NRXN1","ITGB8","CD3G","CD3E","CD79A","PTPRC")
FeaturePlot(int.seu,features = stringr::str_to_title(unlist(features)),order = T,
            reduction = "umap",raster = F)&scale_colour_viridis_c(option = "B")
ggsave("mouse_featureplot2.png",height = 12,width = 16)

DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 3,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP.png",dpi = 300,width = 8,height = 6)
DimPlot(int.seu, group.by = c("technique","treat"),reduction = "umap", label = T,repel = T,label.size = 2,label.box = T,cols = sample(MetBrewer::met.brewer("Klimt",7),2),label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP_subtype.png",dpi = 300,width = 8,height = 6)

# cell.prop <- as.data.frame(prop.table(table(int.seu$celltype)))
# colnames(cell.prop)<- c("sample","proportion")
# label_1 <- paste0(substr(cell.prop$proportion*100,1,4),"%")
# p4=ggplot(cell.prop,aes(sample,proportion,fill=sample))+
#   scale_fill_manual(values = seurat_colors)+
#   geom_bar(stat="identity")+
#   geom_text(aes(label=label_1))+
#   theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
#                    title = element_blank(),
#                    axis.text.x = element_text(angle = 45,hjust = 1))+
#   guides(fill=guide_legend(title=NULL))
# p4

qs::qsave(int.seu,file = "mouse.qs")

#############################
library(stringr)
seurat_colors <- as.character(met.brewer("Klimt", 8))
colors_continuous <- as.character(rev(met.brewer("Hiroshige", 100)))
colors_scale <- as.character(rev(met.brewer("Johnson",100)))
sub_cols <- c("#CF5A79","#544799")

lineage <- read.csv("../NMF/module4.csv",header = T,na.strings = "")
gsc.seu <- subset(int.seu,PiggybacTransposase>0&Akaluc>0&technique=="fresh")

gsc.seu <- AddModuleScore(gsc.seu,features = list(str_to_title(lineage$NPC)),nbin = 30,ctrl = 100,name = "NPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(str_to_title(lineage$AC)),nbin = 30,ctrl = 100,name = "AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(str_to_title(lineage$OPC)),nbin = 30,ctrl = 100,name = "OPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(str_to_title(lineage$cyc)),nbin = 30,ctrl = 100,name = "cycscore")
set.seed(20230829)

lineage.score <- pmax(gsc.seu$NPC1,gsc.seu$AC1)
lineage.score.plot <- lineage.score
lineage.class <- ifelse(gsc.seu$NPC1>gsc.seu$AC1,"NPC","AC")
# lineage.score.plot[lineage.score.plot<0] <- runif(length(lineage.score.plot[lineage.score.plot<0]),0,0.15)
lineage.data <- data.frame(lineage.score,lineage.score.plot,lineage.class)
lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"] <- -lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"]
lineage.data$lineage.class <- ifelse(gsc.seu$OPC1>pmax(gsc.seu$AC1,gsc.seu$NPC1),"OPC",lineage.class)
# lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"]),-0.1,0.1)
# 
# lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"]),-0.1,0.1)
lineage.data$stemness <- gsc.seu$OPC1-lineage.data$lineage.score
lineage.data <- data.frame(lineage.data,gsc.seu@meta.data)
lineage.data$cycling <- ifelse(lineage.data$cycscore1>0,"cycling","non-cycling")

seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
sub_cols <- c("#CF5A79","#544799")

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=lineage.class))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors[c(1,3,5,7,9)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))
ggsave("./2dim-VP-fresh.pdf",height = 3,width = 7.5)

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=cycling))+
  geom_point(size=0.8)+scale_colour_manual(values = c("darkred","#8080804D"))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "top",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )
ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=H3K27M))+
  geom_point(size=0.8)+scale_colour_manual(values = c("darkblue","#8080804D"))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "top",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )
ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=lineage.class))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors[c(1,3,5,7,9)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "top",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )
ggsave("./2dim-VP-fresh.pdf",height = 4,width = 5)


