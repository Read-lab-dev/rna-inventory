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
library(clusterProfiler)
library(enrichplot)
setwd("~/rna/Bulk_Analysis/DIPG")
load("BT17_VP_short_GSEA.Rdata")
DIPG17_GSEA <- gsea_results_df[,-1]
DIPG17_DEG <- resdata
load("BT109_VP_Short_GSEA.Rdata")
BT109_GSEA <- gsea_results_df[,-1]
BT109_DEG <- resdata
overlaps_GSEA <- inner_join(DIPG17_GSEA,BT109_GSEA,by="Description")
DIPG17_DEG <- DIPG17_DEG %>% filter(change=="down")
BT109_DEG <- BT109_DEG %>% filter(change=="down")
overlaps_DEG<- inner_join(DIPG17_DEG[,-1],DIPG17_DEG[,-1],by="gene")
overlaps_DEG <- overlaps_DEG %>% filter(log2FoldChange.x*log2FoldChange.y>0)

##GO analysis
library(org.Hs.eg.db)
overlaps <-bitr(overlaps_DEG$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
enrich.go <- enrichGO(gene = overlaps$ENTREZID,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID',ont = 'ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2,readable = T)
dotplot(enrich.go,showCategory = 30)
ggsave(filename = "DIPG_overlaps_GO_VP_short.pdf",height = 12)
write.csv(enrich.go@result,file = "DIPG_overlaps_GO_VP_short.csv")
write.csv(overlaps_GSEA[,c(1,4,10,13,19)],"DIPG_overlaps_GSEA_VP_short.csv",row.names = F)

enrich.kegg <- enrichKEGG(gene = overlaps$ENTREZID,pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(enrich.kegg,showCategory = 30)
ggsave(filename = "DIPG_overlaps_KEGG_Down_VP_Short.pdf",height = 12)
enrich.kegg <- pairwise_termsim(enrich.kegg)
emapplot(enrich.kegg)
ggsave(filename = "DIPG_overlaps_KEGG_Down_VP_Short_cnet.pdf",width = 12,height = 8)

#######Long Treatment#############
setwd("~/rna/Bulk_Analysis/DIPG")
load("BT17_VP_Long_GSEA.Rdata")
DIPG17_GSEA <- gsea_results_df[,-1]
DIPG17_DEG <- resdata
load("BT109_VP_Long_GSEA.Rdata")
BT109_GSEA <- gsea_results_df[,-1]
BT109_DEG <- resdata
overlaps_GSEA <- inner_join(DIPG17_GSEA,BT109_GSEA,by="Description")

#Up regulated
DIPG17_DEG <- DIPG17_DEG %>% filter(change=="up")
BT109_DEG <- BT109_DEG %>% filter(change=="up")
overlaps_DEG<- inner_join(DIPG17_DEG[,-1],DIPG17_DEG[,-1],by="gene")
overlaps_DEG <- overlaps_DEG %>% filter(log2FoldChange.x*log2FoldChange.y>0)
##GO analysis
library(org.Hs.eg.db)
overlaps <-bitr(overlaps_DEG$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
enrich.go <- enrichGO(gene = overlaps$ENTREZID,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID',ont = 'ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2,readable = T)
dotplot(enrich.go,showCategory = 30)
ggsave(filename = "DIPG_overlaps_GO_UP_VP_Long.pdf",height = 12)
enrich.go <- pairwise_termsim(enrich.go)
emapplot(enrich.go)
ggsave(filename = "DIPG_overlaps_GO_UP_VP_Long_cnet.pdf",width = 12,height = 8)

enrich.kegg <- enrichKEGG(gene = overlaps$ENTREZID,pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(enrich.kegg,showCategory = 30)
ggsave(filename = "DIPG_overlaps_KEGG_UP_VP_Long.pdf",height = 12)
enrich.kegg <- pairwise_termsim(enrich.kegg)
emapplot(enrich.kegg)
ggsave(filename = "DIPG_overlaps_KEGG_UP_VP_Long_cnet.pdf",width = 12,height = 8)


write.csv(enrich.go@result,file = "DIPG_overlaps_GO_UP_VP_Long.csv")
write.csv(overlaps_GSEA[,c(1,4,10,13,19)],"DIPG_overlaps_GSEA_VP_Long.csv",row.names = F)

enrich.go <- setReadable(enrich.go, 'org.Hs.eg.db', 'ENTREZID')
enrich.go <- pairwise_termsim(enrich.go,showCategory = 200)
emapplot(enrich.go,showCategory=30)
ggsave("DIPG_overlaps_GO_UP_VP_Long_cnet.pdf",height = 10,width = 16)

#down regulated
DIPG17_DEG <- DIPG17_DEG %>% filter(change=="down")
BT109_DEG <- BT109_DEG %>% filter(change=="down")
overlaps_DEG<- inner_join(DIPG17_DEG[,-1],DIPG17_DEG[,-1],by="gene")
overlaps_DEG <- overlaps_DEG %>% filter(log2FoldChange.x*log2FoldChange.y>0)
##GO analysis
library(org.Hs.eg.db)
overlaps <-bitr(overlaps_DEG$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
enrich.go <- enrichGO(gene = overlaps$ENTREZID,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID',ont = 'ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2,readable = T)
dotplot(enrich.go,showCategory = 30)
ggsave(filename = "DIPG_overlaps_GO_Down_VP_Long.pdf",height = 12)
write.csv(enrich.go@result,file = "DIPG_overlaps_GO_Down_VP_Long.csv")
enrich.go <- setReadable(enrich.go, 'org.Hs.eg.db', 'ENTREZID')
enrich.go <- pairwise_termsim(enrich.go,showCategory = 200)
emapplot(enrich.go,showCategory=30)
ggsave("DIPG_overlaps_GO_UP_VP_Long_cnet.pdf",height = 10,width = 16)

enrich.kegg <- enrichKEGG(gene = overlaps$ENTREZID,pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(enrich.kegg,showCategory = 30)
ggsave(filename = "DIPG_overlaps_KEGG_DOWN_VP_Long.pdf",height = 12)
enrich.kegg <- pairwise_termsim(enrich.kegg)
emapplot(enrich.kegg)
ggsave(filename = "DIPG_overlaps_KEGG_DOWN_VP_Long_cnet.pdf",width = 12,height = 8)

###########
xx1 <- overlaps$ENTREZID
xx2 <- overlaps$ENTREZID

gcSample <- list(upregulated=xx1,
                 downregualted=xx2)
enrich.go <- compareCluster(gcSample, fun="enrichPathway", OrgDb="org.Hs.eg.db")
enrich.go <- pairwise_termsim(enrich.go)
emapplot(enrich.go)
