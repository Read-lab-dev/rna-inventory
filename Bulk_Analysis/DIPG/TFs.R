rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
setwd("/home/hzg/rna/Bulk_Analysis/DIPG/")
TEAD4 <- read.csv("43301_gene_score_5fold.txt")
TEAD4 <- TEAD4 %>% filter(score>0)
TEAD4 <- TEAD4[!duplicated(TEAD4$gene),]

BT109_VP_s<- read.csv("BT109_VP_short_DEseq2.csv",row.names = 1)
BT109_VP_s <- BT109_VP_s %>% filter(change%in%c("up","down"))

BT109_VP_l<- read.csv("BT109_VP_Long_DEseq2.csv",row.names = 1)
BT109_VP_l <- BT109_VP_l %>% filter(change%in%c("up","down"))

TEAD4_s <-inner_join(TEAD4,BT109_VP_s)
TEAD4_l <-inner_join(TEAD4,BT109_VP_l)

gsea_results_list<- lapply(gsea_results, function(x){
  cat(paste(dim(x@result)),'\n')
  x@result
})
gsea_results_df <- do.call(rbind, gsea_results_list)