###Detect H3K27M transcripts in Bam FIle:
###Zhengang Hu 02-19-2024
rm(list = ls())
setwd("~/rna/BT109VP")
int.seu <-qs::qread("int.seu.qs")

barcode <- read.delim2("H3K27M_allid.txt",header = F)

barcode_WT <- read.delim2("H3K27WT_allid.txt",header = F)

int.seu$barcode <-sapply(strsplit(colnames(int.seu), "_"),tail,1)

H3K27Mcells <- colnames(int.seu)[int.seu$barcode%in%barcode$V1]

H3K27WTcells <- colnames(int.seu)[int.seu$barcode%in%barcode_WT$V1&int.seu$celltype=="DIPG-GSC"]

gfpcells <- colnames(int.seu)[int.seu$gfp=="eGFP+"]
int.seu$gfp <- ifelse(GetAssayData(int.seu)["eGFP",]>0,"eGFP+","eGFP-")

int.seu$H3K27M <- ifelse(int.seu$barcode%in%barcode$V1,"K27M","wild-type")

# int.seu$H3K27WT <- ifelse(int.seu$barcode%in%barcode_WT$V1,"wild-type","undetected")

DimPlot(int.seu, cells.highlight= H3K27Mcells)

qs::qsave(int.seu,file = "int.seu.new.qs")

gsc.seu <-qs::qread("neurosphere.qs")
gsc.seu$barcode <-sapply(strsplit(colnames(gsc.seu), "_"),tail,1)
H3K27Mcells <- colnames(gsc.seu)[gsc.seu$barcode%in%barcode$V1]
DimPlot(gsc.seu, cells.highlight= H3K27Mcells)

table(int.seu$gfp,int.seu$H3K27M)

gsc.seu<- qs::qread("gsc.seu.qs")

barcode_gsc <-sapply(strsplit(colnames(gsc.seu), "_"),tail,1)

H3K27Mcells <- colnames(gsc.seu)[barcode_gsc%in%barcode$V1]
gsc.seu$H3K27M <- ifelse(colnames(gsc.seu)%in%H3K27Mcells,"K27M","wild-type")

DimPlot(gsc.seu, cells.highlight= H3K27Mcells)
