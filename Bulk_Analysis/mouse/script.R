rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(tibble)
library(MetBrewer)
library(ggrepel)
library(ggplot2)
library(ggthemes)
library("RColorBrewer")
library(pheatmap)
library(sva)
library(DESeq2)
###PCA######
library(FactoMineR)
library(factoextra)
gene.pca <- PCA(t(exprSet), ncp = 2, scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(gene.pca,addEllipses = T,ellipse.level=0.95,
             habillage = metadata$Cell,label="ind")
###########DEG##############
load("raw.Rdata")
metadata <- metadata[9:12,]
expr <- exprSet[rowSums(exprSet==0)<0.8*ncol(expr),rownames(metadata)]
dds <- DESeqDataSetFromMatrix(countData=expr, 
                              colData=metadata, 
                              design= ~Cell, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

{res <- results(dds,contrast = c("Cell","PKC","PHC"))
res <- res[order(res$padj),]
resdata <-  as.data.frame(res)
resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 1,
                                 ifelse(resdata$log2FoldChange > 1,"up", "down"), "not"))
resdata$gene <- rownames(resdata)}

write.csv(resdata,"./DEseq2_PKCvsPHC.csv")

#########VolcanoPlot###################
DEG <- na.omit(resdata)
this_tile <- paste0("PKC vs PHC",
                    '\nUp regulated: ',nrow(DEG[DEG$change =='up',]) ,
                    '\nDown regulated: ',nrow(DEG[DEG$change =='down',]))
ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(padj), colour=change)) +
  scale_color_manual(values = c('#99CC00','grey80','#FF9999'))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-1,1),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_label_repel(data = DEG %>% 
                     dplyr::filter(change!="NOT") %>% 
                     top_n(8,log2FoldChange),
                   aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  geom_label_repel(data = DEG %>% 
                     dplyr::filter(change!="NOT") %>% 
                     top_n(8,-log2FoldChange),
                   aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 Fold Change") + ylab("-log10 adjp-value") + ggtitle(this_tile)+
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.position = c(1,0.8),
        legend.justification = c(1, 1))
ggsave(filename = "PKC_VS_PHC_Volcano.png",width  = 6,height = 6,dpi = 300)
###############GSEA#########
{ library(clusterProfiler)
  genelist <- resdata %>% arrange(desc(log2FoldChange))
  gsea_input <- genelist$log2FoldChange
  names(gsea_input) <- rownames(genelist)
  ##Select your Gene Set
  dir="/home/hzg/rna/Bulk_Analysis/MsigDB/mouse/"
  gmts <- list.files(dir,pattern = 'gmt')
  gmts
  #Start GSEA Analysis
  library(GSEABase)
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
}
write.csv(gsea_results_df[,-c(1:2)],file = paste0("GSEA_PKCvsPHC",".csv"))
save(gsea_results,gsea_results_df,resdata,file = paste0("GSEA_PKCvsPHC",".Rdata"))
##########Heatmap#####################
library(tinyarray)
library(ggplot2)
library(cowplot)

test <- DEG %>% filter(change!="NOT") %>% rownames()

fix(submeta)
fac <- factor(submeta$Treat)

draw_heatmap(subexp[test,],fac,show_column_title = T,
             cluster_cols = T,legend =T,
             show_rownames = T,annotation_legend = T,
             n_cutoff = 3,main = celltype,
             color =(rev(met.brewer("Paquin",100,type = "continuous"))))
ggsave(path = paste0("./result/",dataset),filename = paste0("HeatMap_17_",idx,".pdf"),height = 8,width = 7)

########Not Used#####
# sampleDists <- dist(t(assay(vsd)))
# sampleDistMatrix <- as.matrix(sampleDists)
# colnames(sampleDistMatrix) <- NULL
# rownames(sampleDistMatrix) <- paste0(GBM39$Batch,"_",GBM39$Cell)
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,show_rownames = T,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)

# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# df <- as.data.frame(colData(dds)[,c("Cell","Treat")])
# ann_colors <- list(Cell = c("BT17" = "skyblue", "BT109" = "#ff5964", "GBM39" = "#ffe74c"))
# pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=T, annotation_col=df,annotation_colors = ann_colors)