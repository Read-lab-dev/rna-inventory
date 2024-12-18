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
ggsave("Cell.state-by EGFR.pdf")
#######################
cell.prop<-as.data.frame(prop.table(table(hierarchy.result$IGF2BP3,hierarchy.result$cyc)))
colnames(cell.prop)<- c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  scale_fill_manual(values = met.brewer("Klimt",4))+
  geom_bar(stat="identity",position="fill")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
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
p1+p2+p3
ggsave(p1+p2+p3,filename="2dim-expression.png",height = 3,width = 12)

p5=ggplot(data = hierarchy.result %>% arrange(desc(cyc)),aes(x=X,y=Y))+
  geom_point(size=0.5,aes(color=cyc))+
  scale_color_manual(values = c("darkred","#999999"))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))

p6=ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$cyc=="cycling",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_viridis_c(option = "D",name="cycling cell")+scale_alpha_continuous(range = c(0,1),guide=guide_none())
p5+p6
ggsave(p5+p6,filename="2dim-cycling.png",height = 3,width = 9)

p7=ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$RIOK2>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_viridis_c(option = "D",name="RIOK2 density")+scale_alpha_continuous(range = c(0,1),guide=guide_none())
p8=ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$IMP3>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_viridis_c(option = "D",name="IMP3 density")+scale_alpha_continuous(range = c(0,1),guide=guide_none())
p9=ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$MYC>0,],
                 aes(fill=..density..),geom ="raster",contour = F)+scale_fill_viridis_c(option = "D",name="MYC density")+scale_alpha_continuous(range = c(0,1),guide=guide_none())
p7+p8+p9
ggsave(p7+p8+p9,filename="2dim-exp-density.png",height = 3,width = 12)

library(ggnewscale)
ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$cyc=="cycling",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_viridis_c()+scale_alpha_continuous(range = c(0,1),guide=guide_none())+
  new_scale_fill()+
  new_scale("alpha")+
  stat_density2d(data = hierarchy.result[hierarchy.result$RIOK2>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_viridis_c(option = "B",name="RIOK2+")+scale_alpha_continuous(range = c(0,1),guide=guide_none())


ggplot(data=count_dat, mapping=aes(x=cellstate,y=Frac,fill=MYC))+
  geom_bar(stat="identity",width=0.5,position='dodge')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='identity',aes(y=Frac+2,label=paste0(round(Frac,0),"%")),color="black",position = position_dodge(0.5),fontface="bold")+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
ggsave("prop.pdf",width = 5,height = 3)