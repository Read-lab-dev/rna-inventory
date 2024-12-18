library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(gghighlight)
library(MetBrewer)
library(tibble)
Idents(gbm.seu) <- hierarchy.result$cyc
aa <- FindAllMarkers(gbm.seu,test.use = "MAST",logfc.threshold = 0.25)

gbm.seu$cyc <- ifelse(gbm.seu$G1S<0.75&gbm.seu$G2M<0.5,"non-cycling","cycling")

hierarchy.result$RIOK2 <- GetAssayData(gbm.seu)["RIOK2",]
hierarchy.result$IMP3 <-  GetAssayData(gbm.seu)["IGF2BP3",]
hierarchy.result$MYC <-  GetAssayData(gbm.seu)["MYC",]
hierarchy.result$YAP1 <-  GetAssayData(gbm.seu)["YAP1",]
hierarchy.result$TAZ <-  GetAssayData(gbm.seu)["WWTR1",]
hierarchy.result$EGFR <- NULL
pd <- read.csv("neftel_pd.csv",na.strings = "",row.names = 1)[,c(1,5)]
hierarchy.result <- left_join(hierarchy.result,pd,by="patient")
hierarchy.result$patient <- factor(hierarchy.result$patient,levels = pd$patient)
ggplot(data = hierarchy.result,aes(x=X,y=Y,color=cellstate))+
  geom_point(alpha=0.5,size=0.6)+
  scale_color_manual(values = c("#544799","#CF5A79","#D2C564","#5DA373"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  facet_wrap(vars(patient))
ggsave("Cell.state-by EGFR-facet.pdf",height = 6,width = 9)

hierarchy.result <- filter(hierarchy.result,EGFR=="Positive")
p1=ggplot(data = hierarchy.result %>% arrange(RIOK2),aes(x=X,y=Y))+
  geom_point(size=0.5,aes(color=RIOK2))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  scale_color_viridis_c(option = "B")+
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
p2=ggplot(data = hierarchy.result %>% arrange(IMP3),aes(x=X,y=Y))+
  geom_point(size=0.5,aes(color=IMP3))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  scale_color_viridis_c(option = "B")+
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
p3=ggplot(data = hierarchy.result %>% arrange(MYC),aes(x=X,y=Y))+
  geom_point(size=0.5,aes(color=MYC))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  scale_color_viridis_c(option = "B")+
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
p4=ggplot(data = hierarchy.result %>% arrange(YAP1),aes(x=X,y=Y))+
  geom_point(size=0.5,aes(color=YAP1))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  scale_color_viridis_c(option = "B")+
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
p5=ggplot(data = hierarchy.result %>% arrange(TAZ),aes(x=X,y=Y))+
  geom_point(size=0.5,aes(color=TAZ))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  scale_color_viridis_c(option = "B")+
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))

pp1 <- p1|p2|p3|p4|p5
