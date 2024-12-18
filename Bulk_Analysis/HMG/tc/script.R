rm(list = ls())
options(stringsAsFactors = F)
setwd("~/rna/Bulk_Analysis/tc")
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
exprSet <- read.csv(file = './BAM/FeatureCounts.csv',row.names = 1,check.names = F)
####Filter Counts########
table(is.na(exprSet))
colpalette <- sample(rep(met.brewer("Cross",8),20),ncol(exprSet))
boxplot(log2(exprSet+1),col=colpalette)

metadata <- data.frame(row.names = colnames(exprSet),
                       name = rownames(metadata),
                       treat= as.factor(c(rep("ctrl",3),rep("shTX",3))))
save(metadata,exprSet,file = "raw.Rdata")
#################Anno################
anno <- read.delim2("~/rna/human_genome/STAR/geneInfo.tab",skip = 1,header = F)
colnames(anno) <- c("ENSG","gene","type")
anno <- anno[!duplicated(anno$gene),]
anno <- anno %>% filter(type=="protein_coding")
anno <- anno[,1:2]

exprSet <- as.data.frame(exprSet)
exprSet$ENSG <- as.character(rownames(exprSet))
exprSet <- exprSet %>% inner_join(anno) %>% arrange(gene)
rownames(exprSet) <- exprSet$gene
exprSet<- exprSet %>% dplyr::select(-c("ENSG","gene"))
save(exprSet,metadata,file = "rawcounts_coding.Rdata")
###PCA######
library(FactoMineR)
library(factoextra)
gene.pca <- PCA(t(exprSet), ncp = 2, scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(gene.pca,addEllipses = T,ellipse.level=0.95,
             habillage = metadata$treat,label="ind")
###########DEG##############
expr <- exprSet[rowSums(exprSet==0)<0.8*ncol(exprSet),rownames(metadata)]
dds <- DESeqDataSetFromMatrix(countData=expr, 
                              colData=metadata, 
                              design= ~treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

{res <- results(dds,contrast = c("treat","shTX","ctrl"))
  res <- res[order(res$padj),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                   ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
  resdata$gene <- rownames(resdata)
  DEG <- na.omit(resdata)
}

write.csv(resdata,"./DEseq2_shTXvsctrl.csv")

#########VolcanoPlot###################
DEG <- na.omit(resdata)
this_tile <- paste0("shTX vs ctrl",
                    '\nUp regulated: ',nrow(DEG[DEG$change =='up',]) ,
                    '\nDown regulated: ',nrow(DEG[DEG$change =='down',]))
ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(padj), colour=change)) +
  scale_color_manual(values = c('#99CC00','grey80','#FF9999'))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+
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
  theme_clean()+xlim(-3,3)+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.position = c(1,0.8),
        legend.justification = c(1, 1))
ggsave(filename = "shTX_VS_ctrl_Volcano.png",width  = 6,height = 6,dpi = 300)
#############KEGG###########
library(clusterProfiler)
geneList <- DEG %>% filter(change=="down")
geneList <- bitr(gene,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kk_down <- enrichKEGG(gene = geneList$ENTREZID,
                 organism  = 'hsa',
                 pvalueCutoff = 0.05)

###############GSEA#########
{ library(clusterProfiler)
  genelist <- bitr(resdata$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(resdata, by= c("SYMBOL"="gene")) %>% 
    arrange(desc(log2FoldChange))
  gsea_input <- genelist$log2FoldChange
  names(gsea_input) <- genelist$ENTREZID
  ##Select your Gene Set
  dir="/home/hzg/rna/Bulk_Analysis/MsigDB"
  gmts <- list.files(dir,pattern = 'gmt')[-3]
  gmts
  #Start GSEA Analysis
  library(GSEABase)
  gsea_results <- lapply(gmts, function(gmtfile){
    # gmtfile=gmts[2]
    geneset <- read.gmt(file.path(dir,gmtfile)) 
    print(paste0('Now process the ',gmtfile))
    egmt <- GSEA(gsea_input, TERM2GENE=geneset, verbose=FALSE,pvalueCutoff = 0.2)
    head(egmt)
    return(egmt)
  })
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results[[1]] <- setReadable(gsea_results[[1]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results[[2]] <- setReadable(gsea_results[[2]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results[[3]] <- setReadable(gsea_results[[3]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results[[4]] <- setReadable(gsea_results[[4]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results_df <- do.call(rbind, gsea_results_list)
}
write.csv(gsea_results_df[,-c(1:2)],file = paste0("GSEA",".csv"))
save(gsea_results,gsea_results_df,resdata,file = paste0("GSEA",".Rdata"))
enrichplot::gseaplot2(gsea_results[[4]],geneSetID = 1:2,pvalue_table = T)

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