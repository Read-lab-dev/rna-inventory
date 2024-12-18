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
setwd("~/rna/Bulk_Analysis/integrated/202304")
# exprSet <- read.delim2(file = 'FeatureCounts.count',row.names = 1,check.names = F)
# colnames(exprSet) <- gsub("BAM.","A",colnames(exprSet));colnames(exprSet) <- substr(colnames(exprSet),1,13)
# ####Filter Counts########
# table(is.na(exprSet))
# colpalette <- sample(rep(met.brewer("Cross",8),20),ncol(exprSet))
# boxplot(log2(exprSet+0.001),col=colpalette)
# 
# metadata <- read.csv("../BulkSheet.csv",header = T,row.names = 1,na.strings = "")
# metadata <- metadata %>% filter(Cell=="BT17"&Batch=="423")
# 
# exprSet <- exprSet[,rownames(metadata)]
# metadata$Treat<- as.factor(metadata$Treat)
# metadata$Cell<- as.factor(metadata$Cell)
# metadata$Batch <- as.factor(metadata$Batch)
# metadata$Time <- as.factor(metadata$Time)
# save(metadata,exprSet,file = "raw.Rdata")
# #################Anno################
# anno <- read.delim2("~/rna/human_genome/STAR/geneInfo.tab",skip = 1,header = F)
# colnames(anno) <- c("ENSG","gene","type")
# anno <- anno[!duplicated(anno$gene),]
# anno <- anno %>% filter(type=="protein_coding")
# anno <- anno[,1:2]
# 
# exprSet <- as.data.frame(exprSet)
# exprSet$ENSG <- as.character(rownames(exprSet))
# exprSet <- exprSet %>% inner_join(anno) %>% arrange(gene)
# rownames(exprSet) <- exprSet$gene
# exprSet<-exprSet %>% select(-c("ENSG","gene"))
# save(exprSet,metadata,file = "rawcounts_coding.Rdata")
###########################
#######BT17-VP-Paired#######
###########################
load("rawcounts_coding.Rdata")
A20023R-03-01
A20023R-03-02
A20023R-03-03
A20023R-03-04
A20023R-03-05
A20023R-03-06
A20023R-03-07
A20023R-03-08
A20023R-03-09
A20023R-03-19
A20023R-03-20
A20023R-03-21

sample <- clipr::read_clip()
metadata <- metadata[sample,]
metadata$Date <- gsub("/","_",metadata$Date)
metadata$Treat<- factor(metadata$Treat);metadata$Cell<- factor(metadata$Cell);metadata$Batch <- factor(metadata$Batch);metadata$Time <- factor(metadata$Time);metadata$Date <- factor(metadata$Date)
exprSet <- exprSet[,sample]

########PCA#######
dds <- DESeqDataSetFromMatrix(countData=exprSet, 
                              colData=metadata, 
                              design= ~Date+Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Date)
pcaData <- plotPCA(vsd, intgroup=c("Treat","Date","Batch","Full"), returnData=TRUE)
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label="Full",
                  color = "Treat",size = 1, repel = T, 
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95, 
                  main = paste0("Treat"))

###DEG###
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
{res <- results(dds,contrast = c("Treat","sh66","ctrl"),independentFilter=T)
  res <- res[order(res$padj),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                   ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
  resdata$change <- factor(resdata$change,levels = c("up","not","down"))
  resdata$gene <- rownames(resdata)
  table(resdata$change)
  }  

write.csv(resdata,"./DEseq2_BT17_sh66_short.csv")
#########VolcanoPlot###################
gene2show <- c("MAP2","YAP1","WWTR1","DCX","TOP2A","CSPG4","DNMT1","SOX9")

this_tile <- paste0("BT17_sh66_short",
                    '\nUp regulated: ',nrow(na.omit(resdata[resdata$change =='up',])) ,
                    '\nDown regulated: ',nrow(na.omit(resdata[resdata$change =='down',])))
ggplot(data=na.omit(resdata), aes(x=log2FoldChange, y=-log10(pvalue), colour=change)) +
  scale_color_manual(values = c('skyblue','grey80','#FF8889'))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
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
  xlim(c(-3,3))+ylim(c(0,100))+
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.position = c(1, 0.8),
        legend.justification = c(1, 1))
ggsave(filename = "BT17_sh66_Short_Volcano.pdf",width  = 6,height = 4)

###############GSEA#########
resdata <- read.csv("BT17_sh66_DEseq2.csv")
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
  
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
}
write.csv(gsea_results_df,file = paste0("BT17_sh66_GSEA",".csv"))
save(gsea_results,gsea_results_df,file = paste0("BT17_sh66_GSEA",".Rdata"))
#########Heatmap#####################
library(tinyarray)
library(ggplot2)
library(cowplot)

test <- resdata %>% filter(change!="not") %>% rownames()

fac <- factor(metadata$Treat)

draw_heatmap(exprSet[test,],fac,show_column_title = T,
             cluster_cols = T,legend =T,
             show_rownames = F,annotation_legend = T,
             n_cutoff = 3,main = "celltype",
             color =(rev(met.brewer("Paquin",100,type = "continuous"))))
ggsave(filename = paste0("BT17_sh66_Short_HeatMap.pdf"),height = 8,width = 7)

########################
###########################
#######BT17-VP-long-Paired#######
###########################
rm(list = ls())
A20023R-04-05
A20023R-04-06
A20023R-04-11
A20023R-04-13
A20023R-04-15
A20023R-04-16
sample <- clipr::read_clip()
load("~/rna/Bulk_Analysis/integrated/202304/rawcounts_coding.Rdata")
metadata <- metadata[sample,]
metadata$Date <- gsub("/","_",metadata$Date)
metadata$Treat<- factor(metadata$Treat);metadata$Cell<- factor(metadata$Cell);metadata$Batch <- factor(metadata$Batch);metadata$Time <- factor(metadata$Time);metadata$Date <- factor(metadata$Date)
exprSet <- exprSet[,sample]
# metadata$Date <- factor(c(1,1,2,2,3,3))
########PCA#######
dds <- DESeqDataSetFromMatrix(countData=exprSet, 
                              colData=metadata, 
                              design= ~Date+Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Date)
pcaData <- plotPCA(vsd, intgroup=c("Treat","Date","Batch","Full"), returnData=TRUE)
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label="Full",
                  color = "Treat",size = 1, repel = T, 
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95, 
                  main = paste0("Treat"))


######DEG####
dds <- DESeqDataSetFromMatrix(countData=exprSet, 
                              colData=metadata, 
                              design= ~Date+Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
{res <- results(dds,contrast = c("Treat","VP","ctrl"),independentFilter=T)
  res <- res[order(res$padj),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                   ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
  resdata$change <- factor(resdata$change,levels = c("up","not","down"))
  resdata$gene <- rownames(resdata)
  table(resdata$change)}  

write.csv(resdata,"./BT17_VP_Long_DEseq2.csv")
#########VolcanoPlot###################
gene2show <- resdata %>% filter(change!="not")
gene2show <-grep(c("NOTCH.*"),gene2show$gene,value = T)
gene2show <- c("CDH5","OLIG1","OLIG2","TK1","IMP3","DPP4","DBI","GFAP","APOE","DCT")

this_tile <- paste0("BT17_VP_Long",
                    '\nUp regulated: ',nrow(resdata[resdata$change =='up',]) ,
                    '\nDown regulated: ',nrow(resdata[resdata$change =='down',]))
ggplot(data=na.omit(resdata), aes(x=log2FoldChange, y=-log10(pvalue), colour=change)) +
  scale_color_manual(values = c('#FF8889','grey80','skyblue'))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_label_repel(data = resdata[gene2show,],
                   aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 Fold Change") + ylab("-log10 adjp-value") + ggtitle(this_tile)+ylim(0,50)+
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.position = c(1, 0.8),
        legend.justification = c(1, 1))
ggsave(filename = "BT17_VP_Long_Volcano.pdf",width  = 6,height = 4)

###############GSEA#########
resdata <- read.csv("BT17_VP_Long_DEseq2.csv")
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
  
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
}
write.csv(gsea_results_df,file = paste0("BT17_VP_Long_GSEA",".csv"))
save(gsea_results,gsea_results_df,file = paste0("BT17_VP_Long_GSEA",".Rdata"))
#########Heatmap#####################
library(tinyarray)
library(ggplot2)
library(cowplot)

test <- resdata %>% filter(change!="not") %>% rownames()

fac <- factor(metadata$Treat)

draw_heatmap(exprSet[test,],fac,show_column_title = T,
             cluster_cols = T,legend =T,
             show_rownames = F,annotation_legend = T,
             n_cutoff = 3,main = "celltype",
             color =(rev(met.brewer("Paquin",100,type = "continuous"))))
ggsave(filename = paste0("BT17_VP_Long_HeatMap.pdf"),height = 5,width = 3)

########################
###########################
#######BT17-VP-short-Paired#######
###########################
rm(list = ls())
A20023R-03-01
A20023R-03-04
A20023R-03-07
A20023R-03-08
A20023R-03-09
A20023R-03-19
A20023R-03-20
A20023R-03-21
sample <- clipr::read_clip()
load("rawcounts_coding.Rdata")
metadata <- metadata[sample,]
metadata$Date <- gsub("/","_",metadata$Date)
metadata$Treat<- factor(metadata$Treat);metadata$Cell<- factor(metadata$Cell);metadata$Batch <- factor(metadata$Batch);metadata$Time <- factor(metadata$Time);metadata$Date <- factor(metadata$Date)
exprSet <- exprSet[,sample]
# metadata$Date <- factor(c(1,1,2,2,3,3))
########PCA#######
dds <- DESeqDataSetFromMatrix(countData=exprSet, 
                              colData=metadata, 
                              design= ~Date+Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Date)
pcaData <- plotPCA(vsd, intgroup=c("Treat","Date","Batch","Full"), returnData=TRUE)
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label="Full",
                  color = "Treat",size = 1, repel = T, 
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95, 
                  main = paste0("Treat"))


######DEG####
dds <- DESeqDataSetFromMatrix(countData=exprSet, 
                              colData=metadata, 
                              design= ~Date+Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
{res <- results(dds,contrast = c("Treat","VP","ctrl"),independentFilter=T)
  res <- res[order(res$padj),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                   ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
  resdata$change <- factor(resdata$change,levels = c("up","not","down"))
  resdata$gene <- rownames(resdata)
  table(resdata$change)}  

write.csv(resdata,"./BT17_VP_short_DEseq2.csv")
#########VolcanoPlot###################
gene2show <- resdata %>% filter(change!="not")
gene2show <-grep(c("NOTCH.*"),gene2show$gene,value = T)
gene2show <- c("CDH5","OLIG1","OLIG2","TK1","IMP3","DPP4","DBI","GFAP","APOE","DCT")

this_tile <- paste0("BT17_VP_short",
                    '\nUp regulated: ',nrow(resdata[resdata$change =='up',]) ,
                    '\nDown regulated: ',nrow(resdata[resdata$change =='down',]))
ggplot(data=na.omit(resdata), aes(x=log2FoldChange, y=-log10(pvalue), colour=change)) +
  scale_color_manual(values = c('#FF8889','grey80','skyblue'))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_label_repel(data = resdata[gene2show,],
                   aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 Fold Change") + ylab("-log10 adjp-value") + ggtitle(this_tile)+ylim(0,50)+
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.position = c(1, 0.8),
        legend.justification = c(1, 1))
ggsave(filename = "BT17_VP_short_Volcano.pdf",width  = 6,height = 4)

###############GSEA#########
resdata <- read.csv("BT17_VP_short_DEseq2.csv")
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
  
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
}
write.csv(gsea_results_df,file = paste0("BT17_VP_short_GSEA",".csv"))
save(gsea_results,gsea_results_df,file = paste0("BT17_VP_short_GSEA",".Rdata"))
#########Heatmap#####################
library(tinyarray)
library(ggplot2)
library(cowplot)

test <- resdata %>% filter(change!="not") %>% rownames()

fac <- factor(metadata$Treat)

draw_heatmap(exprSet[test,],fac,show_column_title = T,
             cluster_cols = T,legend =T,
             show_rownames = F,annotation_legend = T,
             n_cutoff = 3,main = "celltype",
             color =(rev(met.brewer("Paquin",100,type = "continuous"))))
ggsave(filename = paste0("BT17_VP_short_HeatMap.pdf"),height = 5,width = 3)