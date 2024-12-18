library(ggplot2)
library(ggpubr)
library(ggthemes)
library(patchwork)
library(dplyr)
load("~/rna/sc/IMP3/TCGA.Rdata")
subtype <- read.csv("tcga-pd.csv")
colnames(subtype)[1] <- "PATIENT"
pd <- left_join(pd,subtype,by="PATIENT")
my_comparisons <- list(c("Mesenchymal", "Classical"), 
                        c("Classical", "Proneural"),
                        c("Mesenchymal", "Proneural"))
ggplot(pd,aes(Expression.Subclass,riskscore,fill=Expression.Subclass))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA,alpha=0.9)+
  geom_jitter(width = 0.2,color="black",alpha=0.5)+
  ylab("Expression.Subclass")+
  scale_fill_manual(values = c("#F76D5E", "#FFFFBF", "#72D8FF","purple","red"))+
  stat_compare_means(method="wilcox",comparisons=my_comparisons,
                     label="p.signif")+
  theme_few()+
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text  = element_text(face = "bold"),
        axis.title.x = element_blank())

