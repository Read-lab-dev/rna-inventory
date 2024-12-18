rm(list = ls())
library(clusterProfiler)
library(dplyr)
int.seu <- qs::qread("tumor_adj.qs")

int.seu$cell_orig <- paste0(int.seu$orig.ident,"_",int.seu$celltype)
cellfordeg<-levels(int.seu$celltype)
Idents(int.seu) <- int.seu$celltype.group
DEG_RESULT <- NULL
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(int.seu, ident.1 = paste0(cellfordeg[i],"_CrizS"), 
                         ident.2 = paste0(cellfordeg[i],"_DMSO"), verbose = FALSE,logfc.threshold = 0,test.use = "MAST",max.cells.per.ident = 400)
  CELLDEG$celltype <- cellfordeg[i]
  CELLDEG$symbol <- rownames(CELLDEG)
  DEG_RESULT <- rbind(DEG_RESULT,CELLDEG)
}
write.csv(DEG_RESULT,file = "DEG_result_ALL.csv")
DEG_ALL <-read.csv("D.csv")
###############GSEA#########
for (i in unique(DEG_ALL$celltype)){
  DEG <- DEG_ALL %>% filter(celltype==i)
  ##Creating Gene List for GSEA
  genelist <- bitr(DEG$X,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(DEG, by= c("SYMBOL"="X")) %>% 
    arrange(desc(avg_log2FC))
  gsea_input <- genelist$avg_log2FC
  names(gsea_input) <- genelist$ENTREZID
  ##Select your Gene Set
  dir='/home/hzg/lab/Bulk_analysis/MsigDB/'
  gmts <- list.files(dir,pattern = 'gmt')
  gmts
  #Start GSEA Analysis
  gsea_results <- lapply(gmts, function(gmtfile){
    # gmtfile=gmts[2]
    geneset <- read.gmt(file.path(dir,gmtfile)) 
    print(paste0('Now process the ',gmtfile))
    egmt <- GSEA(gsea_input, TERM2GENE=geneset, verbose=FALSE)
    head(egmt)
    return(egmt)
  })
  gsea_results[[1]] <- setReadable(gsea_results[[1]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results[[2]] <- setReadable(gsea_results[[2]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
  write.csv(gsea_results_df,file = paste0("./GSEA/GSEA_RESULT_",i,".csv") )
  save(gsea_results,gsea_results_df,file = paste0("./GSEA/GSEA_RESULT_",i,".Rdata"))
}
####

inflam <- list(read.table("inflammatory")$V1)
names(inflam)<- "inflammatory"
library(AUCell)
library(clusterProfiler)
cells_rankings <- AUCell_buildRankings(int.seu@assays$RNA@data, nCores=10, plotStats=TRUE)
c5 <- read.gmt("c5.go.bp.v2023.1.Hs.symbols.gmt")
genename <- c("GOBP_ASTROCYTE_DIFFERENTIATION","GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
              "GOBP_OLIGODENDROCYTE_DEVELOPMENT","GOBP_MESENCHYMAL_STEM_CELL_DIFFERENTIATION",
              "GOBP_NEURON_DEVELOPMENT")
geneSets <- lapply(genename, function(x){c5$gene[c5$term == x]})
names(geneSets) <- genename
cells_AUC <- AUCell_calcAUC(inflam, cells_rankings)
set.seed(123)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

DimPlot(int.seu,cells.highlight = cells_assignment$inflammatory$assignment)

int.seu$inflammatory <- as.numeric(getAUC(cells_AUC)["inflammatory", ])
FeaturePlot(int.seu,features = "MET")

int.seu$neuron <- aucs
FeaturePlot(int.seu,features = "neuron")+viridis::scale_color_viridis(option="A")
VlnPlot(int.seu,features = "neuron",group.by = "celltype",split.by = "orig.ident")
