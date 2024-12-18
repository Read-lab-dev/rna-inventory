####It's way much faster performing SCENIC in Python!###
####Do not use R or It only wastes your time!!##
exprMat <- int.seu@assays$RNA@counts
geneskept1 <- ifelse(rowSums(exprMat, na.rm = T)>3*0.01*ncol(exprMat),TRUE,FALSE)
geneskept2 <- ifelse(rowSums(exprMat > 0, na.rm = T)>3*0.01*ncol(exprMat),TRUE,FALSE)
geneskept <- ifelse(geneskept1==TRUE&geneskept2==TRUE,TRUE,FALSE)
exprMat <- exprMat[geneskept,]
data.table::fwrite(as.data.frame(t(as.matrix(exprMat))),file = "/home/hzg/rna/human_genome/scenic_v10/engrafted.csv",row.names = TRUE)

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

inputDir='/home/hzg/rna/human_genome/scenic_v10/neftel/'
inputDir='~/rna/BT109VP/scenic/result_filbin'
scenicLoomPath=file.path(inputDir,'out_SCENIC.loom')
library(SCENIC)
loom <- open_loom(scenicLoomPath) 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
exprMat <- get_dgem(loom)
close_loom(loom)

rownames(regulonAUC)
names(regulons)
int.seu <- qs::qread("~/rna/sc/IMP3/")
int.seu <- subset(int.seu,GBMType=="Pediatric")
int.seu <- subset(int.seu,CellAssignment=="CellAssignment")

int.seu <- qs::qread("~/rna/BT109VP/filbin integration/filbin-tumor.seu.qs")
sub_regulonAUC <- regulonAUC[,colnames(regulonAUC)%in%colnames(int.seu)]
int.seu <- int.seu[,colnames(sub_regulonAUC)]
dim(sub_regulonAUC)
dim(int.seu)
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(int.seu))
int.seu$celltype <- lineage.data$lineage.class
cellClusters <- data.frame(row.names = colnames(int.seu), 
                           cluster = as.character(int.seu$sig))
cellGroup <- data.frame(row.names = colnames(int.seu), 
                        celltype = int.seu$orig.ident)
head(cellClusters)
head(cellGroup)
sub_regulonAUC[1:4,1:4] 

####Motifs
motifEnrichmentFile <- file.path(inputDir, "reg.csv")
motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=1)[-1,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")
colsToShow <- colnames(motifEnrichment)[-c(2,9:11)]
tableSubset <- motifEnrichment[TF=="SOX2"]
viewMotifs(motifEnrichment, colsToShow=colsToShow,
           dbVersion = "v10|10nr",motifCol = "MotifID",regulonTargets =c("highConfAnnot"))
save(motifEnrichment,sub_regulonAUC,regulonAUC,cellGroup,cellClusters,
     file = 'for_rss_and_visual.Rdata')

UMAP = Embeddings(int.seu, reduction = "dim2")
selectedResolution <- "cluster" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(sub_regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- na.omit(t(scale(t(regulonActivity_byCellType), center = T, scale=T)))
# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
library(MetBrewer)
hm <- draw(ComplexHeatmap::Heatmap(rss, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6),col = rev(met.brewer("Hiroshige",100)))) # row font size
hm
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)
viewTable(topRegulators, options = list(pageLength = 10))
metadata <- int.seu@meta.data
# ########AUC###################
# regulonAUC_2<-regulonAUC@assays@data@listData$AUC
# regulons_binaryAUC<-ifelse(regulonAUC_2>0.00849649862615906,1,0)
# regulonAUC@assays@data@listData$AUC <- regulons_binaryAUC
###############Visulization###############
# regulontoplot <- paste0(readClipboard(),"(+)")
# bind_data <- t(scale(assay(regulonAUC[rownames(regulonAUC)%in%regulontoplot,colnames(int.seu)])))
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), cellAnnotation=cellClusters[colnames(sub_regulonAUC), selectedResolution])
## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
plotRSS_oneSet(rss, setName = "RG") # cluster ID

library(SummarizedExperiment)
bind_data <- scale(t(assay(regulonAUC[,colnames(int.seu)])))
bind_data[bind_data>3] <-3
bind_data[bind_data<-3] =-3
int.seu@meta.data <- metadata
int.seu@meta.data = cbind(int.seu@meta.data,bind_data)
# int.seu@meta.data = cbind(int.seu@meta.data,t(assay(sub_regulonAUC["EOMES(+)",colnames(int.seu)])))
# scenic_UMAP <- function(x){
# FeaturePlot(int.seu,features = paste0(x,"(+)"),pt.size = 0.001,
#             raster = T,raster.dpi = c(300,300))+scale_colour_gradientn(colors = rev(met.brewer("Hiroshige",100)))+
#   theme(axis.text = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank())
#   ggsave(filename = paste0("./scenic/result/TF_int_",x,"(+).png"),units = "mm",height = 140,width = 140,dpi = 300)
# }
# scenic_UMAP(x="EOMES")
scenic_UMAP <- function(x){
  FeaturePlot(int.seu,reduction = "dim2",pt.size = 0.2,label.size = 0.5,
              features = paste0(x,"(+)"))&scale_colour_gradientn(colors = rev(met.brewer("Hiroshige",100)))&NoAxes()
}
scenic_UMAP(c("ZEB1","SOX6","SOX2","SOX9","YY1","DDIT3","SOX4","SOX11","PBX1"))
ggsave("scenic_2dim.pdf",height = 14,width = 9)
#######################
gsc.seu <- int.seu
gsc.seu <- subset(int.seu,eGFP>0)
Idents(gsc.seu)<- gsc.seu$orig.ident
gsc.seu <- subset(gsc.seu,idents=c("eng_DMSO","eng_VP_long"))
VlnPlot(gsc.seu,features = rownames(regulonAUC),group.by = "orig.ident",ncol =10)

