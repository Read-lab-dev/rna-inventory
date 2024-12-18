FeatureScatter(gsc.seu,feature1 = "stemness",feature2 = "GSTEM1",shuffle = T)+
  FeatureScatter(gsc.seu,feature1 = "NPC1",feature2 = "GSTEM1",shuffle = T)+
  FeatureScatter(gsc.seu,feature1 = "OPC1",feature2 = "GSTEM1",shuffle = T)+
  FeatureScatter(gsc.seu,feature1 = "AC1",feature2 = "GSTEM1",shuffle = T)


FeaturePlot(gsc.seu,features = "GSTEM1",reduction = "dim2")+scale_color_viridis(option = "d")

FeatureScatter(gsc.seu,feature1 = "stemness",feature2 = "DGC1",shuffle = F)+
  FeatureScatter(gsc.seu,feature1 = "NPC1",feature2 = "DGC1",shuffle = F)+
  FeatureScatter(gsc.seu,feature1 = "OPC1",feature2 = "DGC1",shuffle = F)+
  FeatureScatter(gsc.seu,feature1 = "AC1",feature2 = "DGC1",shuffle = F)

VlnPlot(gsc.seu,features = "stemness",group.by = "celltype",split.by = "orig.ident",flip = T)
