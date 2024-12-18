##Bulk-RNAseq###
options(stringsAsFactors = F)
setwd("./mouse")
library(dplyr)
library(tibble)
library(MetBrewer)
library(ggrepel)
library(ggplot2)
library(ggthemes)
library(patchwork)
library("RColorBrewer")
library(pheatmap)
library(sva)
library(DESeq2)
##
load("rawcounts_coding.Rdata")
###PKC vs PHC##
A20023R-01-09
A20023R-01-10
A20023R-01-11
A20023R-01-12
A20023R-01-13
A20023R-01-14
A20023R-01-15
A20023R-01-16
#####################GBM301############
sample <- clipr::read_clip()
metadata <- metadata[sample,]
metadata$Treat<- factor(metadata$Treat);metadata$Cell<- factor(metadata$Cell)
exprSet <- exprSet[,sample]

########PCA#######
dds <- DESeqDataSetFromMatrix(countData=exprSet, 
                              colData=metadata, 
                              design= ~Cell, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Treat","Cell","Full"), returnData=TRUE)
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label="Full",
                  color = "Cell",size = 1, repel = T, 
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95, 
                  main = paste0("Treat"))

###DEG###
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
{res <- results(dds,contrast = c("Cell","PKC","PHC"),independentFilter=T)
  res <- res[order(res$padj),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                   ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
  resdata$change <- factor(resdata$change,levels = c("up","not","down"))
  resdata$gene <- rownames(resdata)
  table(resdata$change)
}  

write.csv(resdata,"./DEseq2_GBM301_sh75.csv")
#########VolcanoPlot###################
this_tile <- paste0("GBM301_sh75",
                    '\nUp regulated: ',nrow(na.omit(resdata[resdata$change =='up',])) ,
                    '\nDown regulated: ',nrow(na.omit(resdata[resdata$change =='down',])))
ggplot(data=na.omit(resdata), aes(x=log2FoldChange, y=-log10(padj), colour=change)) +
  scale_color_manual(values = c('#FF8889','grey80','skyblue'))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_label_repel(data = resdata %>%
                     dplyr::filter(change=="up") %>%
                     top_n(3,-padj),
                   aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  geom_label_repel(data = resdata %>%
                     dplyr::filter(change=="up") %>%
                     top_n(3,log2FoldChange),
                   aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  geom_label_repel(data = resdata %>%
                     dplyr::filter(change=="down") %>%
                     top_n(3,-log2FoldChange),
                   aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  geom_label_repel(data = resdata %>%
                     dplyr::filter(change=="down") %>%
                     top_n(3,-padj),
                   aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 Fold Change") + ylab("-log10 adjp-value") + ggtitle(this_tile)+
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.position = c(1, 0.8),
        legend.justification = c(1, 1))
ggsave(filename = "GBM301_sh75_Volcano.pdf",width  = 6,height = 4)

###############GSEA#########
# resdata <- na.omit(resdata)

{ library(clusterProfiler)
  genelist <- bitr(resdata$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(resdata, by= c("SYMBOL"="gene")) %>% 
    arrange(desc(log2FoldChange))
  gsea_input <- genelist$log2FoldChange
  names(gsea_input) <- genelist$ENTREZID
  ##Select your Gene Set
  dir="/home/hzg/rna/Bulk_Analysis/MsigDB/"
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
  gsea_results[[1]] <- setReadable(gsea_results[[1]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results[[2]] <- setReadable(gsea_results[[2]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results[[3]] <- setReadable(gsea_results[[3]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
}
write.csv(gsea_results_df[,-c(1:2)],file = paste0("GSEA_GBM301_sh75",".csv"))
save(gsea_results,gsea_results_df,file = paste0("GSEA_GBM301_sh75",".Rdata"))
#########Heatmap#####################
library(tinyarray)
library(ggplot2)
library(cowplot)
fac <- factor(metadata$Treat)
expr_norm <- limma::removeBatchEffect(exprSet, metadata$Date)
test_up <- resdata %>% filter(change!="up") %>% top_n(15,-padj) %>%  rownames()
test_down <- resdata %>% filter(change!="down") %>% top_n(15,-padj) %>%  rownames()
draw_heatmap(expr_norm[c(test_up,test_down),],fac,show_column_title = T,
             cluster_cols = T,legend =T,
             show_rownames = T,annotation_legend = T,
             n_cutoff = 2,main = "celltype",
             color =(rev(met.brewer("Paquin",100,type = "continuous"))))

ggsave(filename = paste0("HeatMap_GBM301_sh75.pdf"),height = 8,width = 7)
expr_norm <- assay(vsd)

save(expr_norm,phenotype,tag,file = "Scissor_GBM301_shRIOK2.Rdata")
