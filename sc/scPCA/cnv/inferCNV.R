rm(list = ls())
gc()
setwd("./cnv")
library(Seurat)
library(MetBrewer)

DimPlot(int.seu)
dat = as.data.frame(GetAssayData(int.seu,slot='counts',assay = 'RNA'))#表达矩阵
dat[1:4,1:4]
dim(dat)

groupinfo = data.frame(v1=colnames(dat),v2=Idents(int.seu))#样本注释文件
head(groupinfo)

library(AnnoProbe)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
head(geneInfor)


dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)

library(data.table)
expFile='expFile.txt'
fwrite(dat,file = expFile,sep = '\t',quote = F,row.names= T)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)

#至此，三个输入文件制作完成，运行infercnv
rm(list=ls())
gc()
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(infercnv)
library(data.table)
library(dplyr)

expFile='expFile.txt' 
groupFiles='groupFiles.txt'  
geneFile='geneFile.txt'

exp <- fread(expFile)
exp <- as.data.frame(exp)
rownames(exp) <- exp$V1
exp <- exp[,-1]
exp <- as.matrix(exp)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exp,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    chr_exclude = c("chrX", "chrY", "chrM"),
                                    ref_group_names=c("Plasma Cell","Cytotoxic T cell"))

infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir='./infercnv_out' ,  #选择输入文件夹
                              cluster_by_groups=T, no_plot=F,
                              cluster_references=F,
                              leiden_resolution=0.02,
                              tumor_subcluster_partition_method="leiden",
                              analysis_mode = "subcluster",
                              hclust_method="ward.D2", plot_steps=F,num_threads = 10 ,denoise=T, HMM=F)
