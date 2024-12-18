###NMF program compairson###
###########Filbin Dataset#############
rm(list = ls())
library(Seurat)
setwd("~/rna/BT109VP/NMF/")
lineage.our <- read.csv("./module4-filbin.csv",header = T,na.strings = "")
lineage.filbin <- read.csv("./lineage2018.csv",header = T,na.strings = "")
lineage.our <- lineage.our[1:30,]
lineage.filbin <-lineage.filbin[1:30,]
gsc.seu <- qs::qread("../filbin/filbin.qs")lineage.our <- read.csv("./module4-filbin.csv",header = T,na.strings = "")

gsc.seu <- subset(gsc.seu,celltype%in%c("K27M","K27M-OC","K27M-GSC","K27M-cGSC"))
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$NPC)),nbin = 30,ctrl = 100,name = "our.NPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$AC)),nbin = 30,ctrl = 100,name = "our.AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$OPC)),nbin = 30,ctrl = 100,name = "our.OPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$cyc)),nbin = 30,ctrl = 100,name = "our.cyc")

gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$OC)),nbin = 30,ctrl = 100,name = "filbin.OC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$AC)),nbin = 30,ctrl = 100,name = "filbin.AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$OPC)),nbin = 30,ctrl = 100,name = "filbin.OPC.shared")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$OPC2)),nbin = 30,ctrl = 100,name = "filbin.OPC.var")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$cyc)),nbin = 30,ctrl = 100,name = "filbin.cyc")

cor.data <- data.frame(gsc.seu@meta.data[,8:16])
colnames(cor.data) <- substr(colnames(cor.data),1,nchar(colnames(cor.data))-1)
library(corrplot)
library(ggplot2)
p.cor <-cor(cor.data)
p.cor <- p.cor[5:9,1:4]
testRes <- cor.mtest(cor.data,conf.level=0.95)
testRes$p <- testRes$p[5:9,1:4]
pdf(file = "corrplot.fibin.pdf", height = 5,width = 5)
corrplot(p.cor,method = 'square',p.mat = testRes$p,
         title = "Filbin.2018.DMG",
         sig.level = c(0.001,0.01,0.05),pch.cex = 0.9,
         insig = "label_sig",tl.col = "black",
         col=rev(COL2('RdYlBu', 100)))
dev.off()

###NMF program compairson###
###########Our Dataset#############
rm(list = ls())
library(Seurat)
setwd("~/rna/BT109VP/NMF/")
lineage.our <- read.csv("./module4.csv",header = T,na.strings = "")
lineage.filbin <- read.csv("./lineage2018_hgnc.csv",header = T,na.strings = "")
lineage.our <- lineage.our[1:30,]
lineage.filbin <-lineage.filbin[1:30,]
gsc.seu <- qs::qread("../gsc.seu.qs")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$NPC)),nbin = 30,ctrl = 100,name = "our.NPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$AC)),nbin = 30,ctrl = 100,name = "our.AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$OPC)),nbin = 30,ctrl = 100,name = "our.OPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$cyc)),nbin = 30,ctrl = 100,name = "our.cyc")

gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$OC)),nbin = 30,ctrl = 100,name = "filbin.OC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$AC)),nbin = 30,ctrl = 100,name = "filbin.AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$OPC)),nbin = 30,ctrl = 100,name = "filbin.OPC.shared")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$OPC2)),nbin = 30,ctrl = 100,name = "filbin.OPC.var")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$cyc)),nbin = 30,ctrl = 100,name = "filbin.cyc")

cor.data <- data.frame(gsc.seu@meta.data[,23:31])
colnames(cor.data) <- substr(colnames(cor.data),1,nchar(colnames(cor.data))-1)
library(corrplot)
library(ggplot2)
p.cor <-cor(cor.data)
p.cor <- p.cor[5:9,1:4]
testRes <- cor.mtest(cor.data,conf.level=0.95)
testRes$p <- testRes$p[5:9,1:4]
pdf(file = "corrplot.our.pdf", height = 5,width = 5)
corrplot(p.cor,method = 'square',p.mat = testRes$p,
         title = "Engrafted Model",
         sig.level = c(0.001,0.01,0.05),pch.cex = 0.9,
         insig = "label_sig",tl.col = "black",
         col=rev(COL2('RdYlBu', 100)))
dev.off()

###########Jessa Dataset#############
rm(list = ls())
library(Seurat)
setwd("~/rna/BT109VP/NMF/")
lineage.our <- read.csv("./module4-filbin.csv",header = T,na.strings = "")
lineage.filbin <- read.csv("./lineage2018.csv",header = T,na.strings = "")
lineage.our <- lineage.our[1:30,]
lineage.filbin <-lineage.filbin[1:30,]
gsc.seu2 <- gsc.seu
gsc.seu <- subset(gsc.seu,Malignant_normal_consensus_Jessa2022=="Malignant")

gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$NPC)),nbin = 30,ctrl = 100,name = "our.NPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$AC)),nbin = 30,ctrl = 100,name = "our.AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$OPC)),nbin = 30,ctrl = 100,name = "our.OPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.our$cyc)),nbin = 30,ctrl = 100,name = "our.cyc")

gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$OC)),nbin = 30,ctrl = 100,name = "filbin.OC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$AC)),nbin = 30,ctrl = 100,name = "filbin.AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$OPC)),nbin = 30,ctrl = 100,name = "filbin.OPC.shared")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$OPC2)),nbin = 30,ctrl = 100,name = "filbin.OPC.var")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage.filbin$cyc)),nbin = 30,ctrl = 100,name = "filbin.cyc")

cor.data <- data.frame(gsc.seu@meta.data[,14:22])
colnames(cor.data) <- substr(colnames(cor.data),1,nchar(colnames(cor.data))-1)
library(corrplot)
library(ggplot2)
p.cor <-cor(cor.data)
p.cor <- p.cor[5:9,1:4]
testRes <- cor.mtest(cor.data,conf.level=0.95)
testRes$p <- testRes$p[5:9,1:4]
pdf(file = "corrplot.jessa.pdf", height = 5,width = 5)
corrplot(p.cor,method = 'square',p.mat = testRes$p,
         title = "Jessa.2022.Pons",
         sig.level = c(0.001,0.01,0.05),pch.cex = 0.9,
         insig = "label_sig",tl.col = "black",
         col=rev(COL2('RdYlBu', 100)))
dev.off()
