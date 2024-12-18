rm(list = ls())
options(stringsAsFactors = F)
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
setwd("~/rna/Bulk_Analysis/oldns")
exprSet <- data.table::fread("GSE122679_raw_counts_GRCh38.p13_NCBI.tsv.gz")


####Filter Counts########
table(is.na(exprSet))
colpalette <- sample(rep(met.brewer("Cross",8),20),ncol(exprSet))
boxplot(log2(exprSet+0.001),col=colpalette)
#################Anno################
anno <- data.table::fread("Human.GRCh38.p13.annot.tsv.gz")[,1:2]
colnames(anno) <- c("GeneID","gene")
anno <- anno[!duplicated(anno$gene),]
exprSet <- as.data.frame(exprSet)
exprSet <- exprSet %>% inner_join(anno) %>% arrange(gene)
rownames(exprSet) <- exprSet$gene
exprSet<-exprSet %>% select(-c("GeneID","gene"))

metadata <- read.table("metadata.txt")
metadata$cell <- as.factor(metadata$cell)
save(exprSet,metadata,file = "rawcounts_coding.Rdata")

########PCA#######
dds <- DESeqDataSetFromMatrix(countData=exprSet, 
                              colData=metadata, 
                              design= ~cell, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Date)
pcaData <- plotPCA(vsd, intgroup=c("cell","id"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label="name",
                      color = "cell",size = 1, repel = T, 
                      ellipse = T, ellipse.type = "norm",
                      ellipse.alpha = 0.1, ellipse.level = 0.95, 
                      main = paste0("Treat")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

###DEG###
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
{res <- results(dds,contrast = c("cell","914","HNPC"),independentFilter=T)
  res <- res[order(res$padj),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 1,
                                   ifelse(resdata$log2FoldChange > 1,"up", "down"), "not"))
  resdata$change <- factor(resdata$change,levels = c("up","not","down"))
  resdata$gene <- rownames(resdata)
  }  

write.csv(resdata,"./DEseq2_914.csv")

###############GSEA#########
{ library(clusterProfiler)
  genelist <- bitr(resdata$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(resdata, by= c("SYMBOL"="gene")) %>% 
    arrange(desc(log2FoldChange))
  gsea_input <- genelist$log2FoldChange
  names(gsea_input) <- genelist$ENTREZID
  dir="/home/hzg/rna/Bulk_Analysis/MsigDB/"
  gmts <- list.files(dir,pattern = 'gmt')[c(1,2,4)]
  gmts
  library(GSEABase)
  gsea_results <- lapply(gmts, function(gmtfile){
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
library(enrichplot)
id <- c("VERHAAK_GLIOBLASTOMA_MESENCHYMAL","VERHAAK_GLIOBLASTOMA_PRONEURAL",
        "VERHAAK_GLIOBLASTOMA_CLASSICAL")
gseaplot2(gsea_results[[1]],geneSetID = id,base_size = 9,
          color = c("#544799","#D16258","#5DA373"),
          rel_heights = c(1.5, 0.5, 0.5),pvalue_table=T)
write.csv(gsea_results_df[,-2],file = paste0("914_GSEA",".csv"),row.names = F)
save(resdata,gsea_results,gsea_results_df,file = paste0("914_GSEA",".Rdata"))

dds <- estimateSizeFactors(dds); 
dat <- as.data.frame(=counts(dds, normalized=TRUE))
data.table::fwrite(dat,file = "expr.csv",row.names = T)

#########Deconvolution########
res <- read.csv("GBMDeconvoluteR.csv")
res$cell <- metadata$cell
res <- reshape2::melt(res)

ggplot(res,aes(cell,value,col=cell))+
  geom_boxplot()+
  facet_wrap(vars(variable))

#######SSGSEA################
library(GSVA) 
geneset <- read.gmt("VERHAAK_GLIOBLASTOMA_MESENCHYMAL.v2023.2.Hs.gmt")
geneset$term <- gsub("VERHAAK_GLIOBLASTOMA_","",geneset$term)
geneset2 <- split(geneset$gene,geneset$term)

gsvaPar <- gsvaParam(dat, geneset2)
gsva.es <- gsva(gsvaPar, verbose=FALSE)
gsva.res <- as.data.frame(t(gsva.es))
gsva.res$cell <- metadata$cell
gsva.res <- reshape2::melt(gsva.res,id="cell")

write.csv(gsva.res,file = "gsva.res.csv")

ggplot(gsva.res,aes(cell,value,fill=cell))+
  geom_boxplot(color="black")+
  facet_wrap(vars(variable))+
  theme_bw()+scale_fill_manual(values = seurat_colors)
