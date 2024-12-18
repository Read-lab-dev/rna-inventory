library(DESeq2)
library(dplyr)
library(ggplot2)
library(tibble)
library("tximport")
library("readr")
library(ggrepel)
library(ggthemes)
####Bulk###
LGG.mtx <- read.table("./raw_LGG/bulk_quant.tsv",header = T)
HGG.mtx <- read.table("./raw_HGG/bulk_quant.tsv",header = T)
exprSet <- cbind(LGG.mtx,HGG.mtx[,-1])
exprSet <- column_to_rownames(exprSet,var = "gene_id")
#################Anno################
ensembl2symbol <- clusterProfiler::bitr(rownames(exprSet),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
ensembl2symbol <- ensembl2symbol[!duplicated(ensembl2symbol$ENSEMBL),]
ensembl2symbol <- ensembl2symbol[!duplicated(ensembl2symbol$SYMBOL),]
exprSet <- exprSet[ensembl2symbol$ENSEMBL,]
identical(rownames(exprSet),ensembl2symbol$ENSEMBL)
rownames(exprSet)<- ensembl2symbol$SYMBOL

metadata <- read.table("bulk-meta-filter.csv",header = T,row.names = 2)
metadata$Treat <- as.factor(metadata$Treat)
exprSet <- exprSet[,rownames(metadata)]
exprSet <- round(exprSet)
save(exprSet,metadata,file = "HGGLGG.Rdata")

dds <- DESeqDataSetFromMatrix(countData=exprSet, 
                              colData=metadata, 
                              design= ~Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Treat"), returnData=TRUE)
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label="name",
                  color = "Treat",size = 1, repel = T, 
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95, 
                  main = paste0("Treat"))
###DEG###
dds <- DESeq(dds)
{res <- results(dds,contrast = c("Treat","HGG","LGG"),independentFilter=T)
  res <- res[order(res$padj),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                   ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
  resdata$change <- factor(resdata$change,levels = c("up","not","down"))
  resdata$gene <- rownames(resdata)
  table(resdata$change)
}  
write.csv(resdata,"./DEseq2_HGGvsLGG.csv")
#########VolcanoPlot###################
gene2show <- c("YAP1","EGFR","TOP2A","TEAD1","TEAD2","CSPG4","HLA-DRA")

this_tile <- paste0("HGG vs LGG",
                    '\nUp regulated: ',nrow(na.omit(resdata[resdata$change =='up',])) ,
                    '\nDown regulated: ',nrow(na.omit(resdata[resdata$change =='down',])))
ggplot(data=na.omit(resdata), aes(x=log2FoldChange, y=-log10(pvalue), colour=change)) +
  scale_color_manual(values = c('skyblue','grey80','#FF8889'))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+  
  geom_label_repel(data = resdata[gene2show,],
                   aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  # geom_label_repel(data = resdata %>% 
  #                    dplyr::filter(change!="NOT") %>% 
  #                    top_n(8,padj),
  #                  aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  # geom_label_repel(data = resdata %>% 
  #                    dplyr::filter(change!="NOT") %>% 
  #                    top_n(8,-padj),
  #                  aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 Fold Change") + ylab("-log10 adjp-value") + ggtitle(this_tile)+
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.position = c(1, 0.8),
        legend.justification = c(1, 1))
ggsave(filename = "./Volcano_HGGvsLGG.csv.pdf",width  = 6,height = 4)

{ library(clusterProfiler)
  genelist <- bitr(resdata$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(resdata, by= c("SYMBOL"="gene")) %>% 
    arrange(desc(log2FoldChange))
  gsea_input <- genelist$log2FoldChange
  names(gsea_input) <- genelist$ENTREZID
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
write.csv(gsea_results_df,file = paste0("HGGvsLGG_GSEA",".csv"))
save(gsea_results,gsea_results_df,resdata,file = paste0("HGGvsLGG_GSEA",".Rdata"))
