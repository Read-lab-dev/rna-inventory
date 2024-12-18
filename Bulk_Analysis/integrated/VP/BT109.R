##BT109

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
BT109 <- metadata %>% filter(Treat%in%c("VP","DMSO","ctrl")&Cell=="BT109")
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
  resLFC$change= as.factor(ifelse(resLFC$pvalue<0.05 & abs(resLFC$log2FoldChange) > 0.5,
                                  ifelse(resLFC$log2FoldChange > 0.5,"UP", "DOWN"), "NOT"))
  resLFC$gene <- rownames(resLFC)
  res <- res[order(rownames(resLFC)),]
  resdata <-  as.data.frame(res)
  resdata$change= as.factor(ifelse(resdata$pvalue<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                   ifelse(resdata$log2FoldChange > 0.5,"UP", "DOWN"), "NOT"))
  resdata$gene <- rownames(resdata)
}

write.csv(resdata,"./DEG_result/DEseq2_BT109_VP.csv")

#########VolcanoPlot###################
DEG <-resdata
DEG$change <- factor(DEG$change,levels=c("UP","NOT","DOWN"))
this_tile <- paste0("BT109",
                    '\nUp regulated: ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nDown regulated: ',nrow(DEG[DEG$change =='DOWN',]))

label2show <- rbind(DEG %>% dplyr::filter(change=="UP") %>% top_n(5,log2FoldChange),
                    DEG %>% dplyr::filter(change=="DOWN") %>% top_n(5,-log2FoldChange),
                    DEG %>% dplyr::filter(change=="UP") %>% top_n(5,-log10(pvalue)),
                    DEG %>% dplyr::filter(change=="DOWN") %>% top_n(5,-log10(pvalue)))

ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(pvalue), colour=change)) +
  scale_color_manual(values = c('#F39B7F','grey80','#4DBBD5'))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_label_repel(data = label2show,
                   aes(x=log2FoldChange, y=-log10(pvalue), colour=change,label=gene))+
  theme_set(theme_set(theme_bw(base_size=20)))+ 
  xlab("log2 Fold Change") + ylab("-log10 p-value") + ggtitle(this_tile)+
  theme_clean()+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.8,0.8))
ggsave(path = paste0("./DEG_result/"),filename = paste0("Volcano","BT109","VP",".png"),width  = 8,height = 6,dpi = 300)
###############GSEA#########
{
  genelist <- bitr(resdata$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(resdata, by= c("SYMBOL"="gene")) %>% 
    arrange(desc(log2FoldChange))
  
  gsea_input <- genelist$log2FoldChange
  
  names(gsea_input) <- genelist$ENTREZID
  
  ##Select your Gene Set
  
  dir='/home/hzg/rna/Bulk_analysis/MsigDB/'
  
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

enrichplot::gseaplot2(gsea_results[[1]],geneSetID = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_HEDGEHOG_SIGNALING","HALLMARK_MYC_TARGETS_V2"),title = "HALLMARK",color = sample(met.brewer("Klimt",7),3))


enrichplot::gseaplot2(gsea_results[[1]],geneSetID = c("HALLMARK_HEDGEHOG_SIGNALING"),title = "HALLMARK_HEDGEHOG",color = "#F39B7F")

enrichplot::gseaplot2(gsea_results[[2]],geneSetID = c("GOBP_T_CELL_MEDIATED_CYTOTOXICITY"),title = "T_CELL_MEDIATED_CYTOTOXICITY",color = "#F39B7F")
enrichplot::gseaplot2(gsea_results[[2]],geneSetID = c("GOBP_SPROUTING_ANGIOGENESIS"),title = "ANGIOGENESIS",color = "#F39B7F")



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
