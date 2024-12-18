#################
### Figure 2i ###
#################
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gghighlight)
library(MetBrewer)
library(tibble)
gsc.seu <- qs::qread("gsc.seu.qs")
gsc.seu <- subset(gsc.seu,idents=c("Eng_DMSO","Eng_VP_L","Eng_VP_S"))
gsc.seu <- subset(gsc.seu,idents=c("BT109_VP","BT109_DMSO"))
seurat_colors <- as.character(met.brewer("Klimt", 8))
colors_continuous <- as.character(rev(met.brewer("Hiroshige", 100)))
colors_scale <- as.character(rev(met.brewer("Johnson",100)))
sub_cols <- c("#CF5A79","#544799")
DefaultAssay(gsc.seu) <- "RNA"
lineage <- read.csv("./NMF/lineage.csv",header = T,na.strings = "")

gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$AC)),nbin = 30,ctrl = 100,name = "AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$OC)),nbin = 30,ctrl = 100,name = "OC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$stemness)),nbin = 30,ctrl = 100,name = "stemness")
set.seed(20240229)
lineage.score <- pmax(gsc.seu$AC1,gsc.seu$OC1)
lineage.score.plot <- lineage.score
lineage.score.plot[lineage.score.plot<0] <- runif(length(lineage.score.plot[lineage.score.plot<0]),0,0.1)
lineage.class <- ifelse(gsc.seu$AC1>gsc.seu$OC1,"AC","OC")
lineage.data <- data.frame(lineage.score,lineage.score.plot,lineage.class)
lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"] <- -lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"]
lineage.data$stemness <- gsc.seu$stemness1-lineage.data$lineage.score
lineage.data <- data.frame(lineage.data,gsc.seu@meta.data)
lineage.data$eng <- ifelse(lineage.data$orig.ident%in%c("BT109_DMSO","BT109_VP"),"unengrafted","engrafted")
lineage.data$vp <- ifelse(lineage.data$orig.ident%in%c("Eng_VP_S","Eng_VP_L","BT109VP"),"VP_Treated","DMSO")
lineage.data$orig.ident <- factor(lineage.data$orig.ident,levels = c("Eng_DMSO","Eng_VP_S","Eng_VP_L"))

seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
sub_cols <- c("#CF5A79","#544799")

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=eng))+
  geom_point(alpha=0.5,size=1)+scale_colour_manual(values = sub_cols)+
  theme_few()+xlab("AC-like <--------------> OC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) 

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=eng))+
  geom_point(alpha=0.3,size=0.8)+scale_colour_manual(values = sub_cols)+
  theme_few()+xlab("AC-like <--------------> OC-like")+gghighlight(unhighlighted_params =list(colour="grey60",alpha=0.1))+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(eng))


ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=vp))+
  geom_point(alpha=0.5,size=1)+scale_colour_manual(values = sub_cols)+
  theme_few()+xlab("AC-like <--------------> OC-like")+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) 

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=orig.ident))+
  geom_point(alpha=0.5,size=0.8)+scale_colour_manual(values = seurat_colors[c(2,5,7)])+
  theme_few()+xlab("AC-like <--------------> OC-like")+gghighlight(unhighlighted_params =list(colour="grey60",alpha=0.2))+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=Phase))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors[c(2,5,7)])+
  theme_few()+xlab("AC-like <--------------> OC-like")+gghighlight(unhighlighted_params =list(colour="grey60",alpha=0.1))+
  theme(
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))

ggplot(data = lineage.data,aes(x=orig.ident,y=stemness,fill=orig.ident))+
  geom_boxplot()+
  scale_fill_manual(values = seurat_colors[c(2,5,7)])+
  ggpubr::stat_compare_means(comparisons=my_comparisons)+ 
  theme_few()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )
ggsave("./2dim/stem_box.png",height = 368,width = 386,units="px",dpi = 72)

ggplot(data = lineage.data,aes(x=orig.ident,y=AC1,fill=orig.ident))+
  geom_boxplot()+
  scale_fill_manual(values = seurat_colors[c(2,5,7)])+
  ggpubr::stat_compare_means(comparisons=my_comparisons)+ 
  theme_few()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )
ggsave("./2dim/AC_box.png",height = 368,width = 386,units="px",dpi = 72)

ggplot(data = lineage.data,aes(x=orig.ident,y=OC1,fill=orig.ident))+
  geom_boxplot()+
  scale_fill_manual(values = seurat_colors[c(2,5,7)])+
  ggpubr::stat_compare_means(comparisons=my_comparisons)+ 
  theme_few()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )
ggsave("./2dim/OC_box.png",height = 368,width = 386,units="px",dpi = 72)

gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$G1.S)),nbin = 30,ctrl = 100,name = "G1S")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$G2.M)),nbin = 30,ctrl = 100,name = "G2M")
ggplot(gsc.seu@meta.data,aes(x=G1S1,y=G2M1,color=Phase))+
  geom_point()

ggplot(data = gsc.seu@meta.data,aes(x=orig.ident,y=G2M1,fill=orig.ident))+
  geom_boxplot()+
  scale_fill_manual(values = seurat_colors[c(2,5,7)])+
  ggpubr::stat_compare_means(comparisons=my_comparisons)+ 
  theme_few()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=Phase))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = c("grey50","red","red"))+
  theme_few()+xlab("AC-like <--------------> OC-like")+
  theme(
    legend.position = "top",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )

####################Neurosphere##########################
gsc.seu <- qs::qread("gsc.seu.new.qs")
gsc.seu <- qs::qread("neurosphere.qs")
gsc.seu <- merge(gsc.seu,ns.seu)
DefaultAssay(gsc.seu) <- "RNA"

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=orig.ident))+
  geom_point(alpha=0.5,size=0.8)+scale_colour_manual(values = seurat_colors[c(2,3,4,5,7)])+
  theme_few()+xlab("AC-like <--------------> OC-like")+gghighlight(unhighlighted_params =list(colour="grey60",alpha=0.2))+
  theme(
    axis.ticks.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))

ggplot(data = lineage.data,aes(x=orig.ident,y=stemness,fill=orig.ident))+
  geom_boxplot()+
  scale_fill_manual(values = seurat_colors[c(2,3,4,5,7)])+
  ggpubr::stat_compare_means()+ 
  theme_few()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )
ggsave("./2dim/stem_box_ns.png",height = 368,width = 386,units="px",dpi = 72)
ggplot(data = lineage.data,aes(x=orig.ident,y=AC1,fill=orig.ident))+
  geom_boxplot()+
  scale_fill_manual(values = seurat_colors[c(2,5,7)])+
  ggpubr::stat_compare_means()+ 
  theme_few()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )
ggsave("./2dim/AC_box_ns.png",height = 368,width = 386,units="px",dpi = 72)

ggplot(data = lineage.data,aes(x=orig.ident,y=OC1,fill=orig.ident))+
  geom_boxplot()+
  scale_fill_manual(values = seurat_colors[c(2,5,7)])+
  ggpubr::stat_compare_means()+ 
  theme_few()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )
ggsave("./2dim/OC_box_ns.png",height = 368,width = 386,units="px",dpi = 72)
