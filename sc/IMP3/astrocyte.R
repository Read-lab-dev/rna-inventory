sig <- read.csv("astrocyte.csv",na.strings = "")
pd <- read.csv("neftel_pd.csv",na.strings = "",row.names = 1)
library(presto)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggthemes)
library(Seurat)
library(clusterProfiler)
library(GSVA)
col.scale <- rev(MetBrewer::met.brewer("Hiroshige",10))
sig$id <- rownames(sig)
geneset <- reshape2::melt(sig,id.vars="id")
geneset <- na.omit(geneset)
geneset$term <- as.character(geneset$variable)
geneset = geneset %>% split(f = .$term, x = .$value)#list

########PseudoBulk##########
expr <- AverageExpression(gbm.seu, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #过滤细胞表达量全为零的基因
expr <- as.matrix(expr)

gs.score.raw <- gsva(expr, gset.idx.list = geneset, 
                     kcdf="Gaussian",
                     method = "ssgsea",
                     parallel.sz=10L)
data <- as.data.frame(t(gs.score.raw))
data$RIOK2 <- ifelse(t(expr)[,"RIOK2"]>median(t(expr)[,"RIOK2"]),"RIOK2-High","RIOK2-Low")
data$IMP3 <- ifelse(t(expr)[,"IGF2BP3"]>median(t(expr)[,"IGF2BP3"]),"IMP3-High","IMP3-Low")
data <- cbind(data[rownames(pd),],pd)[,-8]
data <- reshape2::melt(data)
ggplot(data,aes(variable,value,fill=RIOK2))+
  stat_boxplot(geom = "errorbar",width=0.5)+
  geom_boxplot(outlier.shape = NA,alpha=0.9,width=0.5)+
  geom_jitter(width = 0.2,alpha=0.5,aes(col=RIOK2))+
  ylab("RIOK2")+
  scale_fill_manual(values = c("#F76D5E", "#72D8FF"))+
  stat_compare_means(method="wilcox",label="p.signif")+
  theme_bw()+
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text  = element_text(face = "bold"),
        axis.title.x = element_blank())
ggsave("Astrocye-dif-Neftel-RIOK2.pdf")
