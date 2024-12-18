##########Deconvolution###########
rm(list = ls())
library(Seurat)
library(ggplot2)
setwd("~/rna/BT109VP")
gsc.seu <- qs::qread("gsc.seu.qs")
dat <- as.data.frame(GetAssayData(gsc.seu,slot = "count",assay = "RNA"))
pheno <-data.frame(row.names = rownames(gsc.seu@meta.data),
                   type=gsc.seu$celltype)
dat <- as.data.frame(pheno)
library(data.table)
write.table(dat,file = "./cibersort/count.txt",sep = '\t',quote = F,row.names= T)
write.table(pheno,file = "./cibersort/pheno.txt",sep = '\t',quote = F,row.names= T)

############Cibersortx#####################
library(dplyr)
load("/home/hzg/rna/Bulk_Analysis/integrated/202304/raw.Rdata")
A20023R-03-10
A20023R-03-11
A20023R-03-12
A20023R-03-22
A20023R-03-23
A20023R-03-24
A20023R-04-28
A20023R-04-29
A20023R-04-30
A20023R-04-31
sample <- clipr::read_clip()
metadata <- metadata[sample,]
anno <- read.delim2("~/rna/human_genome/STAR/geneInfo.tab",skip = 1,header = F)
colnames(anno) <- c("ENSG","gene","type")
anno <- anno[!duplicated(anno$gene),]
# anno <- anno %>% filter(type=="protein_coding")
anno <- anno[,1:2]
exprSet <- as.data.frame(exprSet)
exprSet$ENSG <- as.character(rownames(exprSet))
exprSet <- exprSet %>% inner_join(anno) %>% arrange(gene)
rownames(exprSet) <- exprSet$gene
exprSet<-exprSet %>% select(-c("ENSG","gene"))
expr <- exprSet[,rownames(metadata)]
table(is.na(expr))
write.table(expr,file = "./cibersort/BT109_MIXTURE.txt",sep = '\t',quote = F,row.names= T)
######
decov <- read.csv("./cibersort/CIBERSORTx_Results.csv")
metadata <- cbind(metadata,decov)
dat <- metadata
dat$Treat <- gsub("ctrl","DMSO",dat$Treat)
ggplot(dat,aes(Treat,OPC))+geom_boxplot()
wilcox.test(OPC~Treat,data = dat)
