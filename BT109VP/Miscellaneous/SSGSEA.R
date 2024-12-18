###GSEA
library(presto)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggthemes)
library(Seurat)
library(clusterProfiler)
library(GSEABase)
dir='/home/hzg/rna/Bulk_Analysis/MsigDB/'
gmts <- list.files(dir,pattern = 'gmt')
gsc.seu <- qs::qread("../gsc.seu.qs")
seurat.color <- c("#c93f55","#469d76","#3c4b99","#924099","#f7b900","#df9ed4")

aucgene<- wilcoxauc(gsc.seu, 'cluster')
head(aucgene)

#########################################
###############GSEA#########
#Start GSEA Analysis
for (i in levels(gsc.seu)) { 
  gsea_input<- aucgene %>%
    dplyr::filter(group == i) %>%
    arrange(desc(auc)) %>%
    dplyr::select(feature, auc) %>% 
    deframe()
  gsea_results <- lapply(gmts, function(gmtfile){
    # gmtfile=gmts[2]
    geneset <- read.gmt(file.path(dir,gmtfile)) 
    print(paste0('Now process the ',gmtfile))
    egmt <- GSEA(gsea_input, TERM2GENE=geneset, verbose=FALSE)
    head(egmt)
    return(egmt)
  })
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
  save(gsea_results,gsea_results_df,file = paste0("GSEA_RESULT_",i,".Rdata"))
}


#############GSVA################
library(GSVA)
library(msigdbr)
expr <- AverageExpression(gsc.seu, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #过滤细胞表达量全为零的基因
expr <- as.matrix(expr)
geneset <- read.gmt(file.path(dir,gmts[1])) 
geneset$term <- as.character(geneset$term)
geneset = geneset %>% split(f = .$term, x = .$gene)#list
gsva <- gsva(expr, gset.idx.list = geneset, 
             kcdf="Gaussian",
             method = "gsva",
             parallel.sz=10L)

library(pheatmap)

gs <- c("ZHONG_PFC_C1_OPC","ZHONG_PFC_C8_ORG_PROLIFERATING",
        "DESCARTES_FETAL_CEREBRUM_OLIGODENDROCYTES","DESCARTES_FETAL_CEREBRUM_ASTROCYTES")
gs <- c("REACTOME_SIGNALING_BY_WNT","REACTOME_G2_M_DNA_REPLICATION_CHECKPOINT","BIOCARTA_ERBB4_PATHWAY","GOBP_HIPPO_SIGNALING")
pheatmap(gsva, show_colnames = T,angle_col = "45",scale = "none",
         cluster_row = T,cluster_col = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
gsva.data <- as.data.frame(gsva)
gsva.data$name <- rownames(gsva.data)

###SSGSEA###########
col.scale <- rev(MetBrewer::met.brewer("Hiroshige",10))
expr <- GetAssayData(gsc.seu,slot = "data")
expr <- expr[rowSums(expr)>0,]  #过滤细胞表达量全为零的基因
expr <- as.matrix(expr)
geneset <- read.gmt(file.path(dir,gmts[2])) %>% filter(term%in%gs)
geneset$term <- as.character(geneset$term)
geneset = geneset %>% split(f = .$term, x = .$gene)#list

gs.score.raw <- gsva(expr, gset.idx.list = geneset, 
                     kcdf="Gaussian",
                     method = "ssgsea",
                     parallel.sz=10L)
gs.score <- data.frame(t(gs.score.raw)) %>% scale(center = F,scale = T)
metadata <- gsc.seu@meta.data
gsc.seu@meta.data <- metadata
gsc.seu@meta.data <- cbind(gsc.seu@meta.data,gs.score)
gsc.seu$ZHONG_PFC_C3_ASTROCYTE

FeaturePlot(gsc.seu,features = gs,order = T)&scale_colour_gradientn(colors = col.scale)&NoAxes()&theme(title = element_text(size = 8))
VlnPlot(gsc.seu,features = gs[-2])&scale_colour_gradientn(colors = col.scale)&theme(title = element_text(size = 8))
ggsave("./gsea/ssgsea_sig_feature_plot.pdf",height = 4.5,width = 6)
qs::qsave(gsc.seu,file = "../gsc.seu.qs")

FeaturePlot(gsc.seu,pt.size = 0.1,features = c("EGFR","NRG3","TP53","MKI67"),order=T)&NoAxes()&theme(title = element_text(size = 8))
ggsave("Feature_plot1.pdf",height = 4.5,width = 6)

FeaturePlot(gsc.seu,features = gs)&scale_colour_gradientn(colors = col.scale)&NoAxes()&theme(title = element_text(size = 8))
ggsave("ssgsea_sig_GSEA_plot.pdf",height = 4.5,width = 6)
