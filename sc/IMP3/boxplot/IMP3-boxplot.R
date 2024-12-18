library(ggplot2)
library(ggpubr)
library(ggthemes)
library(patchwork)
library(dplyr)
dat <- read.csv("IGF2BP3-CGGA.csv")
dat2 <- read.csv("IGF2BP3-CGGA-His.csv")
my_comparisons <- list( c("Mesenchymal", "Classical"), 
                        c("Classical", "Proneural"),
                        c("Mesenchymal", "Proneural") )
p1 <- ggplot(dat,aes(Subtype,mRNA,fill=Subtype))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA,alpha=0.9)+
  geom_jitter(width = 0.2,color="black",alpha=0.5)+
  ylab("IMP3")+
  scale_fill_manual(values = c("#F76D5E", "#FFFFBF", "#72D8FF"))+
  stat_compare_means(method="wilcox",comparisons=my_comparisons,
                     label="p.signif")+
  theme_few()+
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text  = element_text(face = "bold"),
        axis.title.x = element_blank())

dat2$Histology <- as.factor(ifelse(dat2$Histology=="GBM","GBM","LGG"))
dat2$Histology <- factor(dat2$Histology,levels = c("LGG","GBM"))
p2 <-ggplot(dat2,aes(Histology,mRNA,fill=Histology))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA,alpha=0.9)+
  geom_jitter(width = 0.2,color="black",alpha=0.5)+
  ylab("IMP3")+
  scale_fill_manual(values = c("#72D8FF","#F76D5E"))+
  stat_compare_means(method="wilcox",comparisons =list(c("GBM","LGG")),label="p.signif")+
  theme_few()+
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text  = element_text(face = "bold"),
        axis.title.x = element_blank())
p1+p2
ggsave(p1+p2,file="CGGA-Boxplot.pdf",width = 7,height = 4)

###########rembrandt############
dat <- read.csv("IGF2BP3-REMBRANDT.csv")
dat2 <- read.csv("IGF2BP3-REMBRANDT-His.csv")
my_comparisons <- list( c("Mesenchymal", "Classical"), 
                        c("Classical", "Proneural"),
                        c("Mesenchymal", "Proneural") )
p1 <- ggplot(dat,aes(Subtype,mRNA,fill=Subtype))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA,alpha=0.9)+
  geom_jitter(width = 0.2,color="black",alpha=0.5)+
  ylab("IMP3")+
  scale_fill_manual(values = c("#F76D5E", "#FFFFBF", "#72D8FF"))+
  stat_compare_means(method="wilcox",comparisons=my_comparisons,
                     label="p.signif")+
  theme_few()+
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text  = element_text(face = "bold"),
        axis.title.x = element_blank())

dat2 <- dat2 %>% filter(Histology%in%c("GBM","Astrocytoma","Oligodendroglioma","Mixed glioma"))
my_comparisons= list(c("GBM","Astrocytoma"),
                    c("GBM","Oligodendroglioma"),
                    c("GBM","Mixed glioma"))
dat2$Histology <- factor(dat2$Histology,levels = c("Astrocytoma","Oligodendroglioma","Mixed glioma","GBM"))

p2 <-ggplot(dat2,aes(Histology,mRNA,fill=Histology))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA,alpha=0.9)+
  geom_jitter(width = 0.2,color="black",alpha=0.5)+
  ylab("IMP3")+
  scale_fill_manual(values = c("#72D8FF","#FFFFBF","gray","#F76D5E"))+
  stat_compare_means(method="wilcox",comparisons =my_comparisons,label="p.signif")+
  theme_few()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15,vjust = 1),
        axis.title = element_text(face = "bold"),
        axis.text  = element_text(face = "bold"),
        axis.title.x = element_blank())
p1+p2
ggsave(p1+p2,file="REMBRANDT-Boxplot.pdf",width = 7,height = 4)

#########TCGA#########
dat1 <- read.csv("IGF2BP3-TCGA-NEW.csv")
# dat <- left_join(dat0,dat1,by="Sample")
# dat[is.na(dat$Subtype.y),4] <- dat[is.na(dat$Subtype.y),2]
# dat$Subtype.y <- gsub("Classic-like","Classical",dat$Subtype.y)
# dat$Subtype.y <- gsub("Mesenchymal-like","Mesenchymal",dat$Subtype.y)
# dat <- dat %>% filter(Subtype.y%in%c("Mesenchymal","Classical","Proneural"))
# colnames(dat)[c(3,4)] <- c("mRNA","Subtype")

dat2 <- read.csv("IGF2BP3-TCGA-His.csv")
my_comparisons <- list( c("Mesenchymal", "Classical"), 
                        c("Classical", "Proneural"),
                        c("Mesenchymal", "Proneural") )
my_comparisons <- list( c("ME", "CL"), c("NE", "CL"),
                        c("ME", "NE"),c("PN", "CL"),
                        c("NE", "PN"), c("PN", "ME"))
p1 <- ggplot(dat1,aes(Subtype,mRNA,fill=Subtype))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA,alpha=0.9)+
  geom_jitter(width = 0.2,color="black",alpha=0.5)+
  ylab("IMP3")+
  scale_fill_manual(values = c("#F76D5E", "#FFFFBF", "#72D8FF","green"))+
  stat_compare_means(method="wilcox",comparisons=my_comparisons,
                     label="p.signif")+
  theme_few()+
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text  = element_text(face = "bold"),
        axis.title.x = element_blank())
p1

my_comparisons= list(c("GBM","Astrocytoma"),
                     c("GBM","Oligodendroglioma"),
                     c("GBM","Oligoastrocytoma"))
dat2$Histology <- factor(dat2$Histology,levels = c("Astrocytoma","Oligodendroglioma","Oligoastrocytoma","GBM"))

p2 <-ggplot(dat2,aes(Histology,mRNA,fill=Histology))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA,alpha=0.9)+
  geom_jitter(width = 0.2,color="black",alpha=0.5)+
  ylab("IMP3")+
  scale_fill_manual(values = c("#72D8FF","#FFFFBF","gray","#F76D5E"))+
  stat_compare_means(method="wilcox",comparisons =my_comparisons,label="p.signif")+
  theme_few()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15,vjust = 1),
        axis.title = element_text(face = "bold"),
        axis.text  = element_text(face = "bold"),
        axis.title.x = element_blank())
p1+p2
ggsave(p1+p2,file="IMP3-TCGA-Boxplot.pdf",width = 7,height = 4)

write.csv(dat1,file = "IMP3-TCGA-BOX-value-left-new.csv",row.names = F)
write.csv(dat2,file = "IMP3-TCGA-BOX-value-right.csv",row.names = F)
