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
tmp <- read.delim2(file = 'C:/Users/read_lab/OneDrive - Emory University/bulk/01bulk/data/FeatureCounts.txt',row.names = 1)
exprSet <- tmp

tmp <- read.delim2(file = 'C:/Users/read_lab/OneDrive - Emory University/bulk/02bulk/data/FeatureCounts.txt',row.names = 1)
exprSet <- cbind(exprSet,tmp)

tmp <- read.delim2(file = 'C:/Users/read_lab/OneDrive - Emory University/bulk/03bulk/data/FeatureCounts.txt',row.names = 1)
exprSet <- cbind(exprSet,tmp)

tmp <- read.delim2(file = 'C:/Users/read_lab/OneDrive - Emory University/bulk/04bulk/data/FeatureCounts.txt',row.names = 1)
exprSet <- cbind(exprSet,tmp)

####Filter Counts########
table(is.na(exprSet))

colpalette <- sample(rep(met.brewer("Cross",8),20),ncol(exprSet))
boxplot(log2(exprSet+0.001),col=colpalette)

metadata <- read.csv("BulkSheet.csv",header = T,row.names = 1,na.strings = "")
colnames(exprSet) <- rownames(metadata)
metadata$Treat<- as.factor(metadata$Treat)
metadata$Cell<- as.factor(metadata$Cell)
metadata$Batch <- as.factor(metadata$Batch)
metadata$Time <- as.factor(metadata$Time)
save(metadata,exprSet,file = "raw.Rdata")
#################Anno################
anno <- read.delim2("/home/hzg/rna/Bulk_Analysis/integrated/geneInfo.tab")
anno <- anno[!duplicated(anno$gene),]
# anno <- anno %>% filter(type=="protein_coding")
anno <- anno[,1:2]

exprSet <- as.data.frame(exprSet)
exprSet$ENSG <- as.character(rownames(exprSet))
exprSet <- exprSet %>% inner_join(anno) %>% arrange(gene)
rownames(exprSet) <- exprSet$gene
exprSet<-exprSet %>% dplyr::select(-c("ENSG","gene"))
save(exprSet,metadata,file = "rawcounts_all.Rdata")
###Split####
GBM39 <- metadata %>% filter(Cell=="GBM39")
GBM39[,c("Cell",'Treat','Time',"Date","Batch")] <- lapply(GBM39[,c("Cell",'Treat','Time',"Date","Batch")], factor)
GBM301 <- metadata %>% filter(Cell=="GBM301")
GBM301[,c("Cell",'Treat','Time',"Date","Batch")] <- lapply(GBM301[,c("Cell",'Treat','Time',"Date","Batch")], factor)
BT109 <- metadata %>% filter(Cell=="BT109")
BT109[,c("Cell",'Treat','Time',"Date","Batch")] <- lapply(BT109[,c("Cell",'Treat','Time',"Date","Batch")], factor)
BT17 <- metadata %>% filter(Cell=="BT17")
BT17[,c("Cell",'Treat','Time',"Date","Batch")] <- lapply(BT17[,c("Cell",'Treat','Time',"Date","Batch")], factor)

##Remove Batch Effect####
Remove_batch_effect <- function(celltype){
cellname <- unique(celltype$Cell)
dds <- DESeqDataSetFromMatrix(countData=exprSet[,rownames(celltype)], 
                              colData=celltype, 
                              design= ~ Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
pcaData_before <- plotPCA(vsd, intgroup=c("Treat","Batch"), returnData=TRUE)
percentVar_before <- round(100 * attr(pcaData_before, "percentVar"))

sva_result <- sva::ComBat_seq(exprSet[,rownames(celltype)], batch = celltype$Batch, group = celltype$Treat)
dds <- DESeqDataSetFromMatrix(countData=sva_result, 
                              colData=celltype, 
                              design= ~ Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("Treat","Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p1 <-ggpubr::ggscatter(data = pcaData_before, x = "PC1", y = "PC2",
                  color = "Batch",size = 1, repel = T, 
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95, 
                  main = paste0("Batch_Before_",cellname)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
p2 <- ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",
                    color = "Batch",size = 1, repel = T, 
                    ellipse = T, ellipse.type = "norm",
                    ellipse.alpha = 0.1, ellipse.level = 0.95, 
                    main = paste0("Batch_Before_",cellname)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
p3 <- ggpubr::ggscatter(data = pcaData_before, x = "PC1", y = "PC2",label = "Treat",font.label=c(8,"plain"),
                  color = "Treat",size = 1, repel = T,
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95, 
                  main = paste0("Batch_Before_",cellname)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
p4 <- ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label = 'Treat',font.label=c(8,"plain"),
                  color = "Treat",size = 1, repel = T, 
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95, 
                  main = paste0("Batch_Before_",cellname)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

p1+p2+p3+p4+plot_layout(widths = c(1,1))

ggsave(paste0("./Batch/PCA_",cellname,".pdf"),height = 10,width = 10)
}
Remove_batch_effect(celltype = GBM39)

expr_GBM39 <- sva::ComBat_seq(exprSet[,rownames(GBM39)], batch = GBM39$Batch, group = GBM39$Treat)
expr_GBM301 <- sva::ComBat_seq(exprSet[,rownames(GBM301)], batch = GBM301$Batch, group = GBM301$Treat)
expr_BT109 <- sva::ComBat_seq(exprSet[,rownames(BT109)], batch = BT109$Batch, group = BT109$Treat)
expr_BT17 <- sva::ComBat_seq(exprSet[,rownames(BT17)], batch = BT17$Batch, group = BT17$Treat)
save(expr_GBM39,expr_GBM301,expr_BT109,expr_BT17,BT109,BT17,GBM301,GBM39,file = "RawCounts.Rdata")

sva_expr <- cbind(expr_GBM39,expr_GBM301,expr_BT17,expr_BT109)
sva_expr <- sva_expr[,rownames(metadata)]
dds <- DESeqDataSetFromMatrix(countData=sva_expr, 
                              colData=metadata, 
                              design= ~ Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Cell","Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",
                  color = "Cell",size = 1, repel = T, 
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95, 
                  main = "PCA") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  ggsave("./Batch/PCA.pdf")
########Normalized########
dds <- DESeqDataSetFromMatrix(countData=expr_BT109, 
                              colData=BT109, 
                              design= ~ Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
norm_BT109 <- assay(vsd)
save(norm_GBM39,norm_GBM301,norm_BT109,norm_BT17,BT109,BT17,GBM301,GBM39,file = "NormCounts.Rdata")

###########DEG##############
load("RawCounts.Rdata")
expr <-expr_BT17
expr <- expr[rowSums(expr==0)<0.7*ncol(expr),]
dds <- DESeqDataSetFromMatrix(countData=expr, 
                              colData=BT17, 
                              design= ~ Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

{res <- results(dds,contrast = c("Treat","VP","DMSO"))
res <- res[order(res$padj),]
resdata <-  as.data.frame(res)
resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 1,
                                 ifelse(resdata$log2FoldChange > 1,"up", "down"), "not"))
resdata$gene <- rownames(resdata)}

write.csv(resdata,"./DEG_result/DEseq2_BT109_VP.csv")

#######################################
geneinfo <- read.delim("./geneinfo.tab",header = T)
lnc_RNA <- geneinfo %>% dplyr::filter(type=="lncRNA") %>% select(gene)
miRNA <- geneinfo %>% dplyr::filter(type=="miRNA") %>% select(gene)
mRNA <- geneinfo %>% dplyr::filter(type=="protein_coding") %>% select(gene)

#########VolcanoPlot###################
this_tile <- paste0(idx,
                    '\nUp regulated: ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nDown regulated: ',nrow(DEG[DEG$change =='DOWN',]))
ggplot(data=DEG, aes(x=logFC, y=-log10(adj.P.Val), colour=change)) +
  scale_color_manual(values = c('#FF9999','grey80','#99CC00'))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-1,1),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_label_repel(data = DEG %>% 
                     dplyr::filter(change!="NOT") %>% 
                     top_n(8,logFC),
                   aes(x=logFC, y=-log10(adj.P.Val), colour=change,label=gene))+
  geom_label_repel(data = DEG %>% 
                     dplyr::filter(change!="NOT") %>% 
                     top_n(8,-logFC),
                   aes(x=logFC, y=-log10(adj.P.Val), colour=change,label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 Fold Change") + ylab("-log10 adjp-value") + ggtitle(this_tile)+
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.position = c(1, 0.2),
        legend.justification = c(1, 1))
ggsave(path = paste0("./result/",dataset),filename = paste0(celltype,"_",idx,"_Volcano",".png"),width  = 6,height = 6,dpi = 300)
###############GSEA#########
{
  library(clusterProfiler)
  genelist <- bitr(DEG$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(DEG, by= c("SYMBOL"="gene")) %>% 
    arrange(desc(logFC))
  
  gsea_input <- genelist$logFC
  
  names(gsea_input) <- genelist$ENTREZID
  
  ##Select your Gene Set
  
  dir='C:/Users/read_lab/OneDrive - Emory University/ScRNA_Org/MsigDB'
  
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
write.csv(gsea_results_df,file = paste0("./result/",dataset,"/",celltype,"_",idx,"_GSEA_RESULT",".csv"))
save(gsea_results,file = paste0("./result/",dataset,"/",celltype,"_",idx,"_GSEA_RESULT",".Rdata"))
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