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
setwd("~/rna/BT109VP/scenic/result_jessa/")
inputDir = getwd()
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
close_loom(loom)

int.seu <- qs::qread("../../jessa-tumor.seu.qs")
sub_regulonAUC <- regulonAUC[,intersect(colnames(int.seu),colnames(regulonAUC))]
int.seu <- int.seu[,intersect(colnames(int.seu),colnames(regulonAUC))]
dim(sub_regulonAUC)
dim(int.seu)
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(int.seu))
cellClusters <- data.frame(row.names = colnames(int.seu), 
                           cluster = Idents(int.seu),
                           group = as.character(int.seu$orig.ident))

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
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType[c("TEAD4(+)","STAT1(+)"),]
), center = T, scale=T))
# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)

viewTable(topRegulators, options = list(pageLength = 10))

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
plotRSS_oneSet(rss, setName = "RG") # cluster ID


#########UMAP###########
library(SummarizedExperiment)
library(RColorBrewer)
cols <-  c("grey70","#FDE0DD","#FCC5C0","#FA9FB5","#F768A1","#DD3497","#AE017E","#7A0177","#49006A")
int.seu <- qs::qread("gsc.seu.qs")
bind_data <- as.data.frame(t(getAUC(regulonAUC[,colnames(int.seu)])))
int.seu@meta.data = cbind(int.seu@meta.data,bind_data)
# int.seu@meta.data = cbind(int.seu@meta.data,t(assay(sub_regulonAUC["EOMES(+)",colnames(int.seu)])))
scenic_UMAP <- function(x){
  FeaturePlot(int.seu,features = paste0(x,"(+)"))+scale_colour_gradientn(colors = cols)+
    theme(axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  ggsave(filename = paste0("./TF_int_",x,"(+).png"),units = "mm",height = 140,width = 200,dpi = 300)
}
scenic_UMAP(x="TEAD4")


