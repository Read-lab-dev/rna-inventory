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
##########GBM39shYAP###########
A20023R-03-35
A20023R-03-36
A20023R-03-37
A20023R-03-38
A20023R-04-19
A20023R-04-20
A20023R-04-21
A20023R-04-22
A20023R-04-23
A20023R-04-24
sample <- clipr::read_clip()
exprSet.sub <- exprSet[,sample]
metadata.sub <- metadata[sample,]
metadata.sub$Treat <- factor(metadata.sub$Treat)
dds <- DESeqDataSetFromMatrix(countData=exprSet.sub, 
                              colData=metadata.sub, 
                              design= ~ Date+Treat, tidy = F)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Date","Treat"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label = "group",
                  color = "Treat",size = 1, repel = T, 
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

{res <- results(dds,contrast = c("Treat","sh","GFP"))
  res <- res[order(res$padj),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                   ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
  resdata$gene <- rownames(resdata)
}
DEG <- na.omit(resdata)
gene_name <- c("MYC","YAP1","WWTR1","EGFR","IGFBP3","MET")
ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(padj), colour=change)) +
  scale_color_manual(values = c("blue", "grey", "red"))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_text_repel(data = DEG %>% filter(gene%in%gene_name),aes(x=log2FoldChange, y=-log10(padj),label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 Fold Change") + ylab("-log10 adjp-value") + 
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.position = "top",
        legend.title = element_blank(),
        legend.background=element_blank())
library(clusterProfiler)
downregulated <- DEG %>% filter(change=="down")
eg2np <- bitr(downregulated$gene,fromType='SYMBOL',toType='ENTREZID', OrgDb='org.Hs.eg.db')
kk_down <- enrichKEGG(gene = eg2np$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)
GO_down <- enrichGO(gene = eg2np$ENTREZID,
                      ont = "BP",OrgDb = "org.Hs.eg.db",
                      pvalueCutoff = 0.05)
upregulated <- DEG %>% filter(change=="up")
eg2np <- bitr(upregulated$gene,fromType='SYMBOL',toType='ENTREZID', OrgDb='org.Hs.eg.db')
kk_up <- enrichKEGG(gene = eg2np$ENTREZID,
                      organism = 'hsa',
                      pvalueCutoff = 0.05)
GO_up <- enrichGO(gene = eg2np$ENTREZID,
                    ont = "BP",OrgDb = "org.Hs.eg.db",
                    pvalueCutoff = 0.05)

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

##########Combined shYAP###########
A20023R-03-35
A20023R-03-36
A20023R-03-37
A20023R-03-38
A20023R-04-19
A20023R-04-20
A20023R-04-21
A20023R-04-22
A20023R-04-23
A20023R-04-24
A20023R-03-31
A20023R-03-32
A20023R-03-33
A20023R-03-34

sample <- clipr::read_clip()
exprSet.sub <- exprSet[,sample]
metadata.sub <- metadata[sample,]
metadata.sub$Treat <- factor(metadata.sub$Treat)
dds <- DESeqDataSetFromMatrix(countData=exprSet.sub, 
                              colData=metadata.sub, 
                              design= ~ Date+Treat, tidy = F)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Date","Treat"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label = "group",
                  color = "Treat",size = 1, repel = T, 
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

{res <- results(dds,contrast = c("Treat","sh","GFP"))
  res <- res[order(res$padj),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                   ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
  resdata$gene <- rownames(resdata)
}
DEG <- na.omit(resdata)
gene_name <- c("MYC","YAP1","WWTR1","EGFR","IGFBP3","MET")
ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(padj), colour=change)) +
  scale_color_manual(values = c("blue", "grey", "red"))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_text_repel(data = DEG %>% filter(gene%in%gene_name),aes(x=log2FoldChange, y=-log10(padj),label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 Fold Change") + ylab("-log10 adjp-value") + 
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.position = "top",
        legend.title = element_blank(),
        legend.background=element_blank())
library(clusterProfiler)
downregulated <- DEG %>% filter(change=="down")
eg2np <- bitr(downregulated$gene,fromType='SYMBOL',toType='ENTREZID', OrgDb='org.Hs.eg.db')
kk_down <- enrichKEGG(gene = eg2np$ENTREZID,
                      organism = 'hsa',
                      pvalueCutoff = 0.05)
GO_down <- enrichGO(gene = eg2np$ENTREZID,
                    ont = "BP",OrgDb = "org.Hs.eg.db",
                    pvalueCutoff = 0.05)
upregulated <- DEG %>% filter(change=="up")
eg2np <- bitr(upregulated$gene,fromType='SYMBOL',toType='ENTREZID', OrgDb='org.Hs.eg.db')
kk_up <- enrichKEGG(gene = eg2np$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
GO_up <- enrichGO(gene = eg2np$ENTREZID,
                  ont = "BP",OrgDb = "org.Hs.eg.db",
                  pvalueCutoff = 0.05)

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

##########Combined VP###########
A20023R-03-16
A20023R-03-17
A20023R-03-18
A20023R-03-28
A20023R-03-29
A20023R-03-30
A20023R-03-13
A20023R-03-14
A20023R-03-15
A20023R-03-25
A20023R-03-26
A20023R-03-27
A20023R-04-32
A20023R-04-33
A20023R-04-34
A20023R-04-35


sample <- clipr::read_clip()
exprSet.sub <- exprSet[,sample]
metadata.sub <- metadata[sample,]
metadata.sub$Time <- factor(metadata.sub$Time)
dds <- DESeqDataSetFromMatrix(countData=exprSet.sub, 
                              colData=metadata.sub, 
                              design= ~ Date+Time, tidy = F)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Date","Time"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label = "group",
                  color = "Time",size = 1, repel = T, 
                  ellipse = T, ellipse.type = "norm",
                  ellipse.alpha = 0.1, ellipse.level = 0.95) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

{res <- results(dds,contrast = c("Time","24h","0h"))
  res <- res[order(res$padj),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                   ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
  resdata$gene <- rownames(resdata)
}
DEG <- na.omit(resdata)
gene_name <- c("MYC","YAP1","WWTR1","EGFR","IGFBP3","MET")
ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(padj), colour=change)) +
  scale_color_manual(values = c("blue", "grey", "red"))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_text_repel(data = DEG %>% filter(gene%in%gene_name),aes(x=log2FoldChange, y=-log10(padj),label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 Fold Change") + ylab("-log10 adjp-value") + 
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.position = "top",
        legend.title = element_blank(),
        legend.background=element_blank())
library(clusterProfiler)
downregulated <- DEG %>% filter(change=="down")
eg2np <- bitr(downregulated$gene,fromType='SYMBOL',toType='ENTREZID', OrgDb='org.Hs.eg.db')
kk_down <- enrichKEGG(gene = eg2np$ENTREZID,
                      organism = 'hsa',
                      pvalueCutoff = 0.05)
GO_down <- enrichGO(gene = eg2np$ENTREZID,
                    ont = "BP",OrgDb = "org.Hs.eg.db",
                    pvalueCutoff = 0.05)
upregulated <- DEG %>% filter(change=="up")
eg2np <- bitr(upregulated$gene,fromType='SYMBOL',toType='ENTREZID', OrgDb='org.Hs.eg.db')
kk_up <- enrichKEGG(gene = eg2np$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
GO_up <- enrichGO(gene = eg2np$ENTREZID,
                  ont = "BP",OrgDb = "org.Hs.eg.db",
                  pvalueCutoff = 0.05)

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

p1=ggplot(gsea_results_list[[4]], aes(x = NES, y = reorder(ID,NES), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "D")+
  theme_bw() +
  labs(x = "Normalized Enrichment Score (NES)", 
       y = "Pathway ID",
       fill = "-log10(adjp)") +
  theme(axis.text.y = element_text(size = 8))+
  ggtitle("GSEA_shYAP/TAZ_vs_shGFP")

p2=ggplot(gsea_results_list[[4]], aes(x = NES, y = reorder(ID,NES), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "D")+
  theme_bw() +
  labs(x = "Normalized Enrichment Score (NES)", 
       y = "Pathway ID",
       fill = "-log10(adjp)") +
  theme(axis.text.y = element_text(size = 8))+
  ggtitle("GSEA_VP6h_vs_Ctrl")

p3=ggplot(gsea_results_list[[4]], aes(x = NES, y = reorder(ID,NES), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "D")+
  theme_bw() +
  labs(x = "Normalized Enrichment Score (NES)", 
       y = "Pathway ID",
       fill = "-log10(adjp)") +
  theme(axis.text.y = element_text(size = 8))+
  ggtitle("GSEA_VP24h_vs_Ctrl")
p3
p1/p2/p3
ggsave(p1/p2/p3,file="GSEA_shYAP_VP.pdf",height = 12,width = 8)
