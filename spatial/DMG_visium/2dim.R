gsc.seu <- qs::qread("./seurat/GBM.merge.qs")
lineage.our <- read.csv("~/rna/BT109VP/NMF/module4.csv",header = T,na.strings = "")
# gsc.seu <- subset(gsc.seu,cyc=="cycling")
for (i in colnames(lineage.our)) {
  gsc.seu <- AddModuleScore(gsc.seu,features = list(lineage.our[,i]),
                            nbin = 30,ctrl = 100,name = i)
}
lineage.score <- pmax(gsc.seu$NPC1,gsc.seu$AC1)
lineage.score.plot <- lineage.score
lineage.class <- ifelse(gsc.seu$NPC1>gsc.seu$AC1,"NPC","AC")
# lineage.score.plot[lineage.score.plot<0] <- runif(length(lineage.score.plot[lineage.score.plot<0]),0,0.15)
lineage.data <- data.frame(lineage.score,lineage.score.plot,lineage.class)
lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"] <- -lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"]
lineage.data$lineage.class <- ifelse(gsc.seu$OPC1>pmax(gsc.seu$AC1,gsc.seu$NPC1),"OPC",lineage.class)
# lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"]),-0.1,0.1)
# 
# lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"]),-0.1,0.1)
lineage.data$stemness <- gsc.seu$OPC1-lineage.data$lineage.score
lineage.data <- data.frame(lineage.data,gsc.seu@meta.data)
lineage.data$cycling <- ifelse(lineage.data$cyc1>0,"cycling","non-cycling")

seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
sub_cols <- c("#CF5A79","#544799")

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=lineage.class))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors[c(1,3,5,7,9)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
ggsave("./2dim-class.pdf",height = 3,width = 7.5)

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color= orig.ident))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = met.brewer("Klimt",6))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  facet_wrap(vars(orig.ident),ncol = 6)
ggsave("./2dim-patient.png",height = 3,width = 15)

library(dplyr)
ggplot(data = lineage.data %>%arrange(desc(cycling)),aes(x=lineage.score.plot,y=stemness,color=cycling))+
  geom_point(alpha=0.5,size=0.8)+scale_colour_manual(values =c("darkred","gray60"))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
ggsave("2dim-cycling.pdf",height = 3,width = 7.5)

gsc.seu$celltype <- lineage.data$lineage.class
gsc.seu$cyc <- lineage.data$cycling
SpatialDimPlot(gsc.seu,group.by = "cyc",ncol = 3)&NoLegend()&scale_fill_manual(values = c("red", "skyblue", "green"))
ggsave("spatialdim-cyc-GBM.png",height = 6,width = 12)
SpatialDimPlot(gsc.seu,group.by = "celltype",ncol = 3)&scale_fill_manual(values = c("red", "skyblue", "green"))
ggsave("spatialdim-celltype-GBM.png",height = 6,width = 12)
gsc.seu$celltype[gsc.seu$cyc=="non-cycling"]=NA
gsc.seu$alltype <- lineage.data$lineage.class

qs::qsave(gsc.seu,file = "./seurat/GBM.merge.qs")
