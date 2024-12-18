library(scAB)
library(Scissor)
library(Seurat)
library(ggthemes)
library(ggplot2)
library(ggnewscale)
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
setwd("~/rna/sc/ISHAN/")
source("../IMP3/bulk/scissor_code.R")
###GBM39-VP-Short##
A20023R-01-25
A20023R-01-26
A20023R-01-27
A20023R-01-28
A20023R-03-13
A20023R-03-14
A20023R-03-15
A20023R-03-25
A20023R-03-26
A20023R-03-27
A20023R-04-32
A20023R-04-34
###GBM301-VP-Short##
A20023R-01-29
A20023R-01-30
A20023R-01-31
A20023R-01-32
A20023R-03-16
A20023R-03-17
A20023R-03-18
A20023R-03-28
A20023R-03-29
A20023R-03-30
###GBM39-VP-Long##
A20023R-04-09
A20023R-04-10
A20023R-04-12
A20023R-04-14

load("/home/hzg/rna/Bulk_Analysis/integrated/rawcounts_coding.Rdata")
sample <- clipr::read_clip()
metadata <- metadata[sample,]
metadata$Treat <- gsub("DMSO","ctrl",metadata$Treat)
metadata$Treat<- factor(metadata$Treat);metadata$Cell<- factor(metadata$Cell);metadata$Batch <- factor(metadata$Batch);metadata$Time <- factor(metadata$Time);metadata$Date <- factor(metadata$Date)
exprSet <- exprSet[,sample]
########PCA#######
dds <- DESeqDataSetFromMatrix(countData=exprSet, 
                              colData=metadata, 
                              design= ~Date+Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Batch)
pcaData <- plotPCA(vsd, intgroup=c("Treat","Date","Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label="name",
                     color = "Treat",size = 1, repel = T, 
                     ellipse = T, ellipse.type = "norm",
                     ellipse.alpha = 0.1, ellipse.level = 0.95, 
                     main = paste0("Treat")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  ggpubr::ggscatter(data = pcaData, x = "PC1", y = "PC2",label="name",
                    color = "Batch",size = 1, repel = T, 
                    ellipse = T, ellipse.type = "norm",
                    ellipse.alpha = 0.1, ellipse.level = 0.95, 
                    main = paste0("Batch")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
###DEG###
dds <- DESeq(dds)
{res <- results(dds,contrast = c("Treat","VP","ctrl"),independentFilter=T)
  res <- res[order(res$padj),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                   ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
  resdata$change <- factor(resdata$change,levels = c("up","not","down"))
  resdata$gene <- rownames(resdata)
  table(resdata$change)
}  
write.csv(resdata,"./DEseq2_GBM301_VP_short.csv")
#########VolcanoPlot###################
gene2show <- resdata %>% filter(change!="not")
gene2show <-grep(c("NOTCH.*"),gene2show$gene,value = T)
gene2show <- c("NOTCH1","NOTCH2","NOTCH3","JAG1","PDGFRA","PROM1","CD44","EGFR","CEBPD","JAK3")

this_tile <- paste0("GBM301_VP_short",
                    '\nUp regulated: ',nrow(resdata %>% filter(change=="up")) ,
                    '\nDown regulated: ',nrow(resdata %>% filter(change=="down")))
ggplot(data=na.omit(resdata), aes(x=log2FoldChange, y=-log10(pvalue), colour=change)) +
  scale_color_manual(values = c('skyblue','grey80','#FF8899'))+
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
  xlim(c(-3,3))+
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.position = c(1, 0.8),
        legend.justification = c(1, 1))
ggsave(filename = "GBM301_VP_Short_Volcano.pdf",width  = 6,height = 4)

###############GSEA#########
setwd("~/rna/sc/ISHAN")
resdata <- read.csv("GBM301_VP_Short_DEseq2.csv")
{ library(clusterProfiler)
  genelist <- bitr(resdata$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(resdata, by= c("SYMBOL"="gene")) %>% 
    arrange(desc(log2FoldChange))
  gsea_input <- genelist$log2FoldChange
  names(gsea_input) <- genelist$ENTREZID
  dir="/home/hzg/rna/Bulk_Analysis/MsigDB/"
  gmts <- list.files(dir,pattern = 'gmt')
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
write.csv(gsea_results_df[,-2],file = paste0("GBM301_VP_Short_GSEA",".csv"),row.names = F)
save(resdata,gsea_results,gsea_results_df,file = paste0("GBM301_VP_Short_GSEA",".Rdata"))

######scissor#####
source("./bulk/scissor_code.R")
bulk_dataset <- assay(vsd)
phenotype <- ifelse(metadata$Treat=="VP",1,0)
tag <- c('ctrl', 'VP')
sc_dataset <- run_seurat(gbm.seu,verbose = FALSE)
infos4 <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag,alpha = 0.05,
                  cutoff = 0.2,family = "binomial", Save_file = "Scissor_GBM301_VP_S.RData")
Scissor_select <- rep("unselcted", ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos4$Scissor_pos] <- "VP"
Scissor_select[infos4$Scissor_neg] <- "ctrl"
gbm.seu <- AddMetaData(gbm.seu, metadata = Scissor_select, col.name = "scissor")
gbm.seu$cellstate <- hierarchy.result$cellstate

UMAP_celltype <- DimPlot(gbm.seu, reduction ="dim2",label = T,pt.size=0.0001,group.by = "cellstate")+scale_color_manual(values = c("#CF5A79","#544799","#D2C564","#5DA373"))+theme_few()&NoLegend()
UMAP_scissor <- DimPlot(gbm.seu, reduction = 'dim2', group.by = 'scissor', cols = c('grey','darkred','royalblue'), pt.size=0.0001)+theme_few()
UMAP_celltype+UMAP_scissor
ggsave(UMAP_celltype+UMAP_scissor,file="Scissor_GBM301_VP_Short_select.pdf",height = 3.5,width = 8)

hierarchy.result$scissor <-Scissor_select
ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$scissor=="VP",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+ggtitle("GBM301 on Neftel et.al")+scale_fill_viridis_c(option = "B",name="VP")+scale_alpha_continuous(range = c(0,1),guide=guide_none())+
  new_scale_fill()+
  new_scale("alpha")+
  stat_density2d(data = hierarchy.result[hierarchy.result$scissor=="ctrl",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_viridis_c(option = "D",name="ctrl")+scale_alpha_continuous(range = c(0,1),guide=guide_none())

ggsave(file="Scissor_GBM301_shRIOK2_ctrl.pdf",height = 4,width = 6)

######
load("~/rna/sc/IMP3/Netfel.Rdata")
hierarchy.result$TAZ <- GetAssayData(gbm.seu)["WWTR1",]
hierarchy.result$YAP1 <- GetAssayData(gbm.seu)["YAP1",]
hierarchy.result$TEAD1 <- GetAssayData(gbm.seu)["TEAD1",]
hierarchy.result$TEAD2 <- GetAssayData(gbm.seu)["TEAD2",]
hierarchy.result$TEAD3 <- GetAssayData(gbm.seu)["TEAD3",]
hierarchy.result$TEAD4 <- GetAssayData(gbm.seu)["TEAD4",]

ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(size=1,aes(color=TAZ,alpha=TAZ))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+scale_alpha(range=c(0.7,1),guide = 'none')+
  scale_color_viridis_c(option = "B")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
ggsave("TAZ-neftel.png",dpi = 200,height = 3,width = 5)

ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(size=1,aes(color=YAP1,alpha=YAP1))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+scale_alpha(range=c(0.7,1),guide = 'none')+
  scale_color_viridis_c(option = "B")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
ggsave("YAP1-neftel.png",dpi = 200,height = 3,width = 5)

ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(size=1,aes(color=TEAD3,alpha=TEAD4))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+scale_alpha(range=c(0.7,1),guide = 'none')+
  scale_color_viridis_c(option = "B")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
ggsave("TEAD4-neftel.png",dpi = 200,height = 3,width = 5)

