
DimPlot(int.seu,group.by = "CellAssignment")+FeaturePlot(int.seu,features = "TXNDC16",order = T)

FeaturePlot(gbm.seu,features = "TXNDC16",order = T,reduction = "dim2")

hierarchy.result$TXNDC16 <- GetAssayData(gbm.seu)["TXNDC16",]
ggplot(data = hierarchy.result %>% arrange(TXNDC16),aes(x=X,y=Y))+
  geom_point(size=0.5,aes(color=RIOK2))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  scale_color_viridis_c(option = "A")+
  facet_wrap(vars(patient))+
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))

ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$TXNDC16>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+ggtitle("TXNDC16+GBM cell Density on Neftel et.al")+scale_fill_viridis_c(option = "B",name="TXNDC16")+scale_alpha_continuous(range = c(0,1),guide=guide_none())
