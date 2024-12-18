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
setwd("~/rna/Bulk_Analysis/HMG")
exprSet <- read.csv("./FeatureCounts-single.csv",row.names = 1)
####Filter Counts########
table(is.na(exprSet))

colpalette <- sample(rep(met.brewer("Cross",8),20),ncol(exprSet))
metadata <- data.frame(row.names = colnames(exprSet),
                       group = rep(c("Ctrl","CoCulture","TMZ","oeTX"),each=4),
                       name= colnames(exprSet))
save(metadata,exprSet,file = "raw.Rdata")
#################Anno################
anno <- read.delim2("/home/hzg/rna/Bulk_Analysis/integrated/geneInfo.tab")
anno <- anno[!duplicated(anno$gene),]
anno <- anno %>% filter(type=="protein_coding")
anno <- anno[,1:2]
# anno$ENSG <- unlist(stringr::str_split(anno$ENSG,"[.]",simplify=T))[,1]
exprSet <- as.data.frame(exprSet)
exprSet$ENSG <- as.character(rownames(exprSet))
exprSet <- exprSet %>% inner_join(anno) %>% arrange(gene)
rownames(exprSet) <- exprSet$gene
exprSet<-exprSet %>% dplyr::select(-c("ENSG","gene"))
save(exprSet,metadata,file = "rawcounts_all.Rdata")

exprSet <- exprSet[rowSums(exprSet==0)<0.8*ncol(exprSet),]
# gene_detect <- resdata %>% filter(abs(log2FoldChange)>=11) %>% rownames()

exprSet_filter <- exprSet[,-c(1,6,10,14)]
metadata_filter <- metadata[-c(1,6,10,14),]
dds <- DESeqDataSetFromMatrix(countData=exprSet_filter, 
                              colData=metadata_filter, 
                              design= ~ group, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)

# library(PCAtools)
# library(patchwork)
# p <- pca(assay(vsd), metadata = metadata_filter, removeVar = 0.1)
# screeplot(p, axisLabSize = 18, titleLabSize = 22)
# biplot(p,colby = "group")
# plotloadings(p, labSize = 3)
# p$loadings
# dat <- reshape2::melt(exprSet_filter[c("TXNDC5","IFIT3","MX1","TAGLN","IGFBP5","AKR1B1","AIF1"),])
# dat <- reshape2::melt(exprSet_filter[c("CD79B","ARG1","CXCL13","PDGFB","VEGFA","TPSAB1"),])
# dat <- reshape2::melt(exprSet_filter[c("HLA-DR","ARG1","CXCL13","PDGFB","VEGFA","TPSAB1"),])
# dat$Var2 <- rep(metadata_filter$group,each=6)
# ggplot(dat, aes(x = as.factor(Var2), y = value)) +
#   geom_boxplot(aes(fill = Var2), position = position_dodge(0.9)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_wrap(vars(gene),scales = "free")


pcaData <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",
                  color = "group",size = 1, repel = T,
                  ellipse = F, ellipse.type = "norm",label="name",
                  ellipse.alpha = 0.1, ellipse.level = 0.7, mean.point = F,
                  main = "PCA") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  ggforce::geom_mark_ellipse(aes(fill = group,
                                 color = group),expand = 0.0001) +
  theme(legend.position = 'top')
ggsave("./PCA.pdf",width = 4,height = 4)
########Normalized########
norm <- assay(vsd)
boxplot(log2(norm),palette=colpalette)
###########DEG##############
dds <- DESeq(dds)
resultsNames(dds) 
{res <- results(dds,contrast = c("group","oeTX","CoCulture"))
# res <- lfcShrink(dds,coef = 4)
resdata <-  as.data.frame(res)
resdata <- resdata[order(resdata$padj),]
resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                 ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
resdata$gene <- rownames(resdata)
DEG <- na.omit(resdata)
}
write.csv(DEG,"./DEseq2_CoCulture_vs_Ctrl.csv")
write.csv(DEG,"./DEseq2_oe_TX_vs_CoCulture.csv")
write.csv(DEG,"./DEseq2_TMZ_vs_CoCulture.csv")
#########VolcanoPlot###################
this_tile <- paste0("oeTXNDC16 vs Ctrl",
                    '\nUp regulated: ',nrow(DEG[DEG$change =='up',]) ,
                    '\nDown regulated: ',nrow(DEG[DEG$change =='down',]))
ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(padj), colour=change)) +
  scale_color_manual(values = c("blue", "grey", "red"))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_text_repel(data = DEG %>% 
                     dplyr::filter(change=="up") %>% 
                     top_n(8,log2FoldChange),
                   aes(x=log2FoldChange, y=-log10(padj),label=gene))+
  geom_text_repel(data = DEG %>% 
                     dplyr::filter(change=="down") %>% 
                     top_n(8,(-log2FoldChange)),
                   aes(x=log2FoldChange, y=-log10(padj),label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 Fold Change") + ylab("-log10 adjp-value") + ggtitle(this_tile)+xlim(c(-2,2))+ylim(c(0,30))+
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.position = "top",
        legend.title = element_blank(),
        legend.background=element_blank())
ggsave(filename = "oeTX_Volcano.png",width= 6,height = 6,dpi = 300)
###############GSEA#########
{library(clusterProfiler)
  genelist <- DEG %>% arrange(desc(log2FoldChange))
  gsea_input <- genelist$log2FoldChange
  names(gsea_input) <- genelist$gene
  dir='/home/hzg/rna/Bulk_Analysis/MsigDB/'
  gmts <- list.files(dir,pattern = 'gmt')
  gmts
  #Start GSEA Analysis
  library(GSEABase)
  gsea_results <- lapply(gmts, function(gmtfile){
    # gmtfile=gmts[2]
    geneset <- read.gmt(file.path(dir,gmtfile)) 
    print(paste0('Now process the ',gmtfile))
    egmt <- GSEA(gsea_input, TERM2GENE=geneset, verbose=FALSE,nPermSimple = 10000)
    head(egmt)
    return(egmt)
  })
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
}
gsea_results_df_Coculture <- gsea_results_df
gsea_results_df_oeTX <- gsea_results_df
gsea_results_df_TMZ <- gsea_results_df

write.csv(gsea_results_df_Coculture[,-c(1:2)],"GSEA_coculture.csv")
write.csv(gsea_results_df_oeTX[,-c(1:2)],"GSEA_oeTX.csv")
write.csv(gsea_results_df_TMZ[,-c(1:2)],"GSEA_TMZ.csv")
save(gsea_results_df_oeTX,gsea_results_list,file = "oeTX-GSEA.Rdata")
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


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- pcaData$name
rownames(sampleDistMatrix) <- pcaData$name
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDistMatrix,show_rownames = T,cluster_rows=T,cluster_cols = T,
         clustering_distance_rows=sampleDists,
        clustering_distance_cols=sampleDists,
        col=colors)

# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# df <- as.data.frame(colData(dds)[,c("Cell","Treat")])
# ann_colors <- list(Cell = c("BT17" = "skyblue", "BT109" = "#ff5964", "GBM39" = "#ffe74c"))
# pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=T, annotation_col=df,annotation_colors = ann_colors)

###HMC3 outside####
expr_hm3 <- data.table::fread("GSE200354_raw_counts_GRCh38.p13_NCBI.tsv.gz")
expr_hm3 <- as.data.frame(expr_hm3)
ref <- data.table::fread("Human.GRCh38.p13.annot.tsv.gz")
expr_hm3$GeneID <- ref$Symbol
metadata2 <- data.frame(row.names = colnames(expr_hm3)[-1],
                        group = rep(c("GCNY","NC","NORMAL"),each=3),
                        name= colnames(expr_hm3)[-1])
exprSet$GeneID <- rownames(exprSet)

exprSet_merge <- inner_join(exprSet,expr_hm3)
meta_merge <- rbind(metadata,metadata2)

exprSet_merge <- tibble::column_to_rownames(exprSet_merge,"GeneID")

exprSet_merge <- exprSet_merge[,-1]
meta_merge <- meta_merge[-1,]
exprSet_merge <- exprSet_merge[rowSums(exprSet_merge==0)<0.8*ncol(exprSet_merge),]
meta_merge$batch <- c(rep("B1",15),rep("B2",9))
# gene_detect <- resdata %>% filter(abs(log2FoldChange)>=11) %>% rownames()
dds <- DESeqDataSetFromMatrix(countData=exprSet_merge, 
                              colData=meta_merge, 
                              design= ~batch+group, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd),batch = meta_merge$batch)
pcaData <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",
                  color = "group",size = 1, repel = T, 
                  ellipse = F, ellipse.type = "norm",label="name",
                  ellipse.alpha = 0.1, ellipse.level = 0.7, mean.point = F,
                  main = "PCA") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  ggforce::geom_mark_ellipse(aes(fill = group,
                                 color = group),expand = 0.0001) +
  theme(legend.position = 'top')
