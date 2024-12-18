rm(list = ls())
gc()
setwd("~/rna/BT109VP/Rcnv")
library(Seurat)
library(MetBrewer)

org.seu <- qs::qread("/home/hzg/rna/BT109VP/int.seu.NEW1.qs")
org.seu <- subset(org.seu,celltype=="pGluN")
org.seu <- subset(org.seu,downsample=400)
int.seu <- subset(int.seu,downsample=400)

int.seu <- merge(int.seu,org.seu)
int.seu <- JoinLayers(int.seu)
dat = as.data.frame(GetAssayData(int.seu,slot='counts',assay = 'RNA'))#表达矩阵
dat[1:4,1:4]
dim(dat)

groupinfo = data.frame(v1=colnames(dat),v2=Idents(int.seu))#样本注释文件
head(groupinfo)

library(AnnoProbe)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
geneInfor$chr <- factor(geneInfor$chr,c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX","chrY","chrM"))
geneInfor=geneInfor[with(geneInfor, order(chr,start)),c(1,4:6)]
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
                                    ref_group_names="pGluN")
options(scipen = 100)
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir='./infercnv' ,  #
                              cluster_by_groups=T, no_plot=F,
                              cluster_references=F,
                              tumor_subcluster_partition_method="leiden",
                              analysis_mode = "subcluster",
                              hclust_method="ward.D2", plot_steps=F,num_threads = 10 ,denoise=T, HMM=T)

infercnv::plot_cnv(infercnv_obj2, #
                   plot_chr_scale = T, #
                   cluster_by_groups=T,cluster_references=F,
                   output_filename = "better_plot",output_format = "png",
                   custom_color_pal=color.palette(c("#2596be","white","#ef382e"), c(100,100))) #保存为pdf文件#改颜色

qs::qsave(infercnv_obj2,file="infercnv_data.qs")
# 
# color.palette(c("#8DD3C7","white","#BC80BD"), c(100,100))
groupings <- read.delim2("./infercnv_out_5/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings")
CNV_status <- read.delim2("./infercnv_out_5/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat")

aa <- read.delim2("./infercnv_out_5/infercnv_subclusters.observation_groupings.txt",sep = " ")
aa$name <- rownames(aa)
aa$name <- gsub("_1","",aa$name)
bb <- left_join(int.seu@meta.data,aa)
View(bb)
int.seu$leiden <- aa$Annotation.Group
DimPlot(int.seu,group.by = "leiden")
table(gsc.seu$celltype,gsc.seu$leiden)

qs::qsave(gsc.seu,"../gsc.seu.qs")

