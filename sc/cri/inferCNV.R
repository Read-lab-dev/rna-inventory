rm(list = ls())
gc()
library(Seurat)
library(MetBrewer)
seurat_colors <- as.character(met.brewer("Klimt", 20))
colors_continuous <- as.character(rev(met.brewer("Hiroshige", 100)))
colors_scale <- as.character(rev(met.brewer("Johnson",100)))
sub.seu <- subset(int.seu,downsample=800)

dat = as.data.frame(sub.seu@assays$RNA$counts)#表达矩阵
dat[1:4,1:4]
dim(dat)

groupinfo = data.frame(v1=colnames(dat),v2=sub.seu$seurat_clusters)#样本注释文件
head(groupinfo)

library(AnnoProbe)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
geneInfor$chr <- factor(geneInfor$chr,c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX","chrY","chrM"))
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
options(scipen = 100)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exp,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    chr_exclude = c("chrX", "chrY", "chrM"),
                                    ref_group_names=NULL)

infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir='./infer_dou' ,  #选择输入文件夹
                              cluster_by_groups=T, no_plot=F,
                              analysis_mode = "subcluster",
                              hclust_method="ward.D2", plot_steps=F,num_threads = 10 ,denoise=T, HMM=F)


infercnv::plot_cnv(infercnv_obj2, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = "plot",output_format = "png",
                   custom_color_pal=color.palette(c("#2596be","white","#ef382e"), c(10,10)))  #保存为pdf文件#改颜色


qs::qsave(infercnv_obj2,file="infercnv_data.qs")
