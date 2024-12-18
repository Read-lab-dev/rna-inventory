astrocyte <- read.csv("../../CopyOfastrocyte.csv",na.strings = "")
astrocyte$id <- rownames(astrocyte)
astrocyte <- na.omit(reshape2::melt(astrocyte,id.vars="id"))[,-1]
colnames(astrocyte) <- c("sig","gene")
load("~/rna/sc/IMP3/bulk/bulk-figure/GSEA_GBM39_shRIOK2.Rdata")
load("~/rna/sc/IMP3/bulk/bulk-figure/GSEA_GBM39_sh77.Rdata")
load("~/rna/sc/IMP3/bulk/bulk-figure/GSEA_GBM301_shRIOK2.Rdata")
load("~/rna/sc/IMP3/bulk/bulk-figure/GSEA_GBM301_sh77.Rdata")
resdata2 <- left_join(resdata,astrocyte)
resdata2 <- resdata2[!is.na(resdata2$sig),] %>% dplyr::group_by(sig)%>% arrange(log2FoldChange,.by_group = T)
resdata2$gene <- factor(resdata2$gene,levels = resdata2$gene)


p1=ggplot(resdata2)+
  geom_bar(stat = "identity",position = "dodge",
           aes(x=gene,y=log2FoldChange,fill=sig))+
  scale_fill_manual(values =  c("#5DA373","#CF5A79","#D2C564"))+
  ggtitle("GBM39 shRIOK2 vs shCtrl")+
  theme_bw()+ylab("log2 Fold Change")+
  theme(legend.position = "top",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5))
p2=ggplot(resdata2)+
  geom_bar(stat = "identity",position = "dodge",
           aes(x=gene,y=log2FoldChange,fill=sig))+
  scale_fill_manual(values =  c("#5DA373","#CF5A79","#D2C564"))+
  ggtitle("GBM39 shIMP3 vs shCtrl")+
  theme_bw()+ylab("log2 Fold Change")+
  theme(legend.position = "top",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5))
p3=ggplot(resdata2)+
  geom_bar(stat = "identity",position = "dodge",
           aes(x=gene,y=log2FoldChange,fill=sig))+
  scale_fill_manual(values =  c("#5DA373","#CF5A79","#D2C564"))+
  ggtitle("GBM301 shRIOK2 vs shCtrl")+
  theme_bw()+ylab("log2 Fold Change")+
  theme(legend.position = "top",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5))
p4=ggplot(resdata2)+
  geom_bar(stat = "identity",position = "dodge",
           aes(x=gene,y=log2FoldChange,fill=sig))+
  scale_fill_manual(values =  c("#5DA373","#CF5A79","#D2C564"))+
  ggtitle("GBM301 shIMP3 vs shCtrl")+
  theme_bw()+ylab("log2 Fold Change")+
  theme(legend.position = "top",axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5))
(p1+p2)/(p3+p4)
ggsave("Bulk-astrocyte-barplot.pdf",width = 12,height = 10)
