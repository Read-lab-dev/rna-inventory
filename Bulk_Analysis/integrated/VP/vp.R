rm(list = ls())
options(stringsAsFactors = F)
library(clusterProfiler)
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
###########################
# res <- results(dds,contrast = c("Batch","323","123"))
# # library("FactoMineR")
# # library("factoextra")
# res.pca <- PCA(t(assay(vsd)), graph = FALSE)
# fviz_pca_ind(res.pca,
#              geom.ind = "point", # show points only (nbut not "text")
#              col.ind = BT109$Treat, # color by groups
#              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#              addEllipses = TRUE, # Concentration ellipses
#              legend.title = "Groups"
# )
# vsd <- vst(dds,blind = T)
# plotPCA(vsd, intgroup=c("Batch"))
# mm <- model.matrix(~Treat, colData(vsd))
# assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=vsd$Batch, design=mm)
# plotPCA(vsd, intgroup=c("Batch"))
#########BT109##########
BT109 <- metadata %>% filter(Treat%in%c("VP","DMSO")&Cell=="BT109")
BT109[,c("Cell",'Treat','Time',"Date","Batch")] <- lapply(BT109[,c("Cell",'Treat','Time',"Date","Batch")], factor)
expr_BT109 <- exprSet[,rownames(BT109)]
dds <- DESeqDataSetFromMatrix(countData=expr_BT109, 
                              colData=BT109, 
                              design= ~ Batch+Treat, tidy = F)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

{
  res <- results(dds,contrast = c("Treat","VP","DMSO"))
  resLFC <- as.data.frame(lfcShrink(dds, contrast = c("Treat","VP", "DMSO"),type="ashr"))
  resLFC$change= as.factor(ifelse(resLFC$padj<0.05 & abs(resLFC$log2FoldChange) > 1,
                                   ifelse(resLFC$log2FoldChange > 1,"UP", "DOWN"), "NOT"))
  resLFC$gene <- rownames(resLFC)
  res <- res[order(rownames(res)),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 1,
                                   ifelse(resdata$log2FoldChange > 1,"UP", "DOWN"), "NOT"))
  resdata$gene <- rownames(resdata)
  }

write.csv(resdata,"./DEG_result/DEseq2_BT109_VP.csv")

#########VolcanoPlot###################
DEG <-resLFC[!is.na(resdata$change),]
this_tile <- paste0("BT109",
                    '\nUp regulated: ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nDown regulated: ',nrow(DEG[DEG$change =='DOWN',]))
ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(pvalue), colour=change)) +
  scale_color_manual(values = c('#99CC00','grey80','#FF9999'))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-1,1),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_label_repel(data = DEG %>% 
                     dplyr::filter(change=="UP") %>% 
                     top_n(5,-padj),
                   aes(x=log2FoldChange, y=-log10(pvalue), colour=change,label=gene))+
  geom_label_repel(data = DEG %>% 
                     dplyr::filter(change=="DOWN") %>% 
                     top_n(5,-padj),
                   aes(x=log2FoldChange, y=-log10(pvalue), colour=change,label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+ 
  xlab("log2 Fold Change") + ylab("-log10 p-value") + ggtitle(this_tile)+
  theme_clean()+
  ylim(c(0,100))+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.position = c(1, 0.8),
        legend.justification = c(1, 1))
ggsave(path = paste0("./DEG_result/"),filename = paste0("Volcano","BT109","VP",".png"),width  = 6,height = 6,dpi = 300)
###############GSEA#########
{
  genelist <- bitr(resLFC$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(resLFC, by= c("SYMBOL"="gene")) %>% 
    arrange(desc(log2FoldChange))
  
  gsea_input <- genelist$log2FoldChange
  
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

write.csv(gsea_results_df,file = "./DEG_result/GSEA_BT109_VP.csv")
save(gsea_results,resdata,file = "./DEG_result/GSEA_BT109_VP.Rdata")
##########Heatmap#####################
library(tinyarray)
library(ggplot2)
library(cowplot)

test <- DEG %>% filter(change=="UP") %>% top_n(5,wt=lfcSE) %>% rownames()
test <- c(test, DEG %>% filter(change=="DOWN") %>% top_n(15,wt=lfcSE) %>% rownames())
fac <- factor(BT109$Treat)
vsd <- vst(dds,blind = T)
mat <- assay(vsd)
mm <- model.matrix(~Treat, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
assay(vsd) <- mat


draw_heatmap(mat[test,],fac,show_column_title = T,
             cluster_cols = T,legend =T,
             show_rownames = T,annotation_legend = T,
             n_cutoff = 3,main = "BT109",
             color =(rev(met.brewer("Paquin",100,type = "continuous"))))
ggsave(filename = "./DEG_result/HeatMap_BT109.pdf",height = 8,width = 7)
