rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
setwd("/home/hzg/rna/BT109VP/scenic/counts")
inputDir = "/home/hzg/rna/BT109VP/scenic/counts"
scenicLoomPath=file.path(inputDir,'out_SCENIC.loom')
library(SCENIC)
loom <- open_loom(scenicLoomPath) 
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
close_loom(loom)

int.seu <- qs::qread("../../int.seu.NEW1.qs")
int.seu <- qs::qread("../../gsc.seu.qs")
# regulonAUC <- read.csv("binarization.csv",row.names = 1)
sub_regulonAUC <- regulonAUC[,intersect(colnames(int.seu),colnames(regulonAUC))]
int.seu <- int.seu[,intersect(colnames(int.seu),colnames(regulonAUC))]
dim(sub_regulonAUC)
dim(int.seu)
Idents(int.seu) <- int.seu$celltype
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(int.seu))
# Idents(int.seu) <- paste0(int.seu$orig.ident,"_",int.seu$celltype)
cellClusters <- data.frame(row.names = colnames(int.seu), 
                           cluster = Idents(int.seu))
                           # group = paste0(int.seu$orig.ident,"_",Idents(int.seu)))

head(cellClusters)
sub_regulonAUC[1:4,1:4] 

selectedResolution <- "cluster" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType),center = T, scale=T))

# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",row_names_gp=grid::gpar(fontsize=8))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)
viewTable(topRegulators, options = list(pageLength = 10))
select_regulon <- topRegulators %>% group_by(CellType) %>% arrange(desc(RelativeActivity)) %>% top_n(15,wt=RelativeActivity) %>% arrange(CellType)
select_regulon <- c("DDIT3(+)", "YY1(+)","ATF4(+)","REST(+)","XBP1(+)","ZNF281(+)","ASCL1(+)","FOXM1(+)","RFXANK(+)","TEAD4(+)","ZBTB18(+)","PBX1(+)","HLTF(+)","SOX11(+)","SOX4(+)","POU3F2(+)","RB1(+)","ZEB1(+)","TCF7L2(+)","SOX6(+)")

regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType),function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[select_regulon,], name="Regulon activity",row_names_gp=grid::gpar(fontsize=6)))

#############Motifs###############
motifEnrichment <- data.table::fread("./reg.csv", header=T, skip=1)[-1,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")
colsToShow <- colnames(motifEnrichment)[-c(2,9:11)]
tableSubset <- motifEnrichment[TF=="SOX6"]
viewMotifs(motifEnrichment, colsToShow=colsToShow,
           dbVersion = "v10|10nr",motifCol = "MotifID",regulonTargets =c("highConfAnnot"))

#############RSS#########
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[colnames(regulonAUC), selectedResolution])
## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
plotRSS_oneSet(rss, setName = c("RG")) # cluster ID

library(ggplot2)
library(ggrepel)
library(dplyr)
plotRSS_all <- function (rss, setName, n = 5)
{
  rssThisType <- melt(rss)
  rssThisType <- rssThisType %>% group_by(Var2) %>% arrange(dplyr::desc(value),.by_group = T)
  rssThisType <- rssThisType %>% group_by(Var2)%>% mutate(rank = row_number(dplyr::desc(value)))
  rssThisType[rssThisType$rank>n,1] <- NA
  ggplot(rssThisType, aes(x = rank, y = value)) +
    geom_point(color = "#758d71", size = 1) + 
    geom_label_repel(aes(label = Var1,colour=factor(Var1)),
                     box.padding = 0.25, point.padding = 0.25, segment.color = "grey50", 
                     na.rm = TRUE) + theme_classic()+
    facet_wrap(vars(Var2),scales="free",nrow = 2)&NoLegend()+ 
    theme(strip.background = element_rect(fill = "#c8b7a6"),
          strip.text = element_text(color = "#675643",face = "bold"))
}
plotRSS_all(rss)
ggsave("RSSplot.all.pdf",height = 8,width = 16)
#########UMAP###########
library(SummarizedExperiment)
library(RColorBrewer)
cols <-  c("grey70","#D8C9C8","#EAD4D2","#FDE0DD","#FCC5C0","#FA9FB5","#F768A1","#DD3497","#AE017E","#7A0177","#49006A")
dat <- scale(t(getAUC(regulonAUC[,colnames(int.seu)])))
bind_data <- as.data.frame(t(getAUC(regulonAUC[,colnames(int.seu)])))
int.seu@meta.data = cbind(int.seu@meta.data,dat)
# int.seu@meta.data = cbind(int.seu@meta.data,t(assay(sub_regulonAUC["EOMES(+)",colnames(int.seu)])))
scenic_UMAP <- function(x){
  FeaturePlot(int.seu,features = paste0(x,"(+)"),order = T)+scale_colour_gradientn(colors = cols,"Regulon")+
    theme(axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
}
scenic_UMAP(x="SOX4")
ggsave("UMAP+SOX4+Regulon.png",height = 4,width = 6)
seurat_colors <- c("#CF5A79","#544799","#D77B5A","#5DA373")

dat <- regulonActivity_byCellType_Scaled[,]
dat_melt <- melt(dat)
dat_melt <- dat_melt %>% group_by(Var2)%>% mutate(rank = row_number(dplyr::desc(value)))
id <- unique(filter(dat_melt,rank<=5) %>% select(Var1))
dat <- regulonActivity_byCellType_Scaled[id$Var1,]
dat_melt <- melt(dat)

library(dendextend)
library(ggplotify)
dat_melt$Var1 <- factor(dat_melt$Var1, levels = id$Var1)
dat_melt$Var2 <- factor(dat_melt$Var2, levels = levels(id$Var2))
ggplot(dat_melt,aes(x=Var1,y=Var2,fill=value))+
  geom_tile(color="white", size=1)+scale_fill_gradient2(low="#003366", high="#990033", mid="white")+
  ggthemes::theme_few()+
  labs(fill='Activity')+
  theme(axis.text.x = element_text(angle=-90,hjust = -0.2,
        colour = rep(seurat_colors,each=5),face = "bold"),
        axis.text.y = element_text(face = "bold",colour = seurat_colors),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank())
ggsave("scenic_heatmap.pdf",height = 4)
