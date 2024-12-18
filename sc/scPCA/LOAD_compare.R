rm(list = ls())
LGG.seu <- qs::qread("scPCA_LGG.qs")
str <- paste0("CU",clipr::read_clip())
LGG.seu <- subset(LGG.seu,orig.ident%in%str)
LGG.seu$class <- "LGG"
HGG.seu <- qs::qread("scPCA_HGG.qs")
HGG.seu$class <- "HGG"
HGG.seu$percent.mt <- HGG.seu$subsets_mito_percent
str <- clipr::read_clip()
HGG.seu <- subset(HGG.seu,orig.ident%in%str)
Idents(HGG.seu) <- HGG.seu$orig.ident
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(HGG.seu)
HGG.seu <- RenameIdents(HGG.seu, new.cluster.ids)
Idents(HGG.seu) <- HGG.seu$celltype
qs::qsave(HGG.seu,file = "HGG.seu")

#######################
int.seu <- merge(LGG.seu,HGG.seu)
int.seu <- SCTransform(int.seu,vars.to.regress = "percent.mt")
int.seu = RunPCA(int.seu,verbose=F)
# int.seu = harmony::RunHarmony(int.seu,group.by.vars="orig.ident",max_iter=20)
int.seu <- RunUMAP(int.seu,reduction = "pca",dims = 1:30)
int.seu = FindNeighbors(int.seu,reduction = "pca",dims=seq(30),verbose=F)
int.seu = FindClusters(int.seu,res=0.5,verbose=F)
DimPlot(int.seu,reduction = "umap",label = T)
seurat_colors <- c("#DF9ED4","#406E89","#D77B5A","#D2C564","#5DA373","#CF5A79","#544799", "#924099")
p1=DimPlot(int.seu,reduction = "umap",label = T, label.box = T,group.by = "celltype",repel = T,label.size = 4,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)&NoLegend()

FeaturePlot(int.seu,features = c("YAP1","WWTR1","TEAD1","TEAD2","TEAD3","TEAD4"),order = T)
p2=FeaturePlot(int.seu,features = c("YAP1","TEAD1"),split.by = "class",order = T)&scale_color_viridis_b(option = "D",alpha = 1)
p1/p2

ggsave(p1/p2,filename="HGGvsLGG_UMAP.png",dpi = 200,height = 12,width = 8)
qs::qsave(int.seu,file = "int.seu.qs")

DimPlot(int.seu,group.by = "tumor")

gsc.seu <- qs::qread("int.seu.qs")
gsc.seu <- subset(gsc.seu,celltype=="Neoplastic")

lineage <- read.csv("../../BT109VP/NMF/module4.csv",header = T,na.strings = "")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$NPC)),nbin = 25,ctrl = 100,name = "NPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$AC)),nbin = 25,ctrl = 100,name = "AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$OPC)),nbin = 25,ctrl = 100,name = "OPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$cyc)),nbin = 25,ctrl = 100,name = "cycscore")
lineage.score <- pmax(gsc.seu$NPC1,gsc.seu$AC1)
lineage.score.plot <- lineage.score
lineage.class <- ifelse(gsc.seu$NPC1>gsc.seu$AC1,"NPC","AC")
lineage.data <- data.frame(lineage.score,lineage.score.plot,lineage.class)
lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"] <- -lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"]
lineage.data$lineage.class <- ifelse(gsc.seu$OPC1>pmax(gsc.seu$AC1,gsc.seu$NPC1),"OPC",lineage.class)
lineage.data$stemness <- gsc.seu$OPC1-lineage.data$lineage.score
lineage.data$cycscore <- gsc.seu$cycscore1
lineage.data$lineage.class <- ifelse(lineage.data$cycscore>0,"cyc",lineage.data$lineage.class)
lineage.data <- data.frame(lineage.data,gsc.seu@meta.data)
seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
lineage.data$YAP1 <- GetAssayData(gsc.seu)["YAP1",]
lineage.data$TEAD1 <- GetAssayData(gsc.seu)["TEAD1",]
ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=lineage.class))+
  geom_point(alpha=1,size=0.1)+scale_colour_manual(values = seurat_colors[c(1,3,5,7,9)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(class))
ggsave("2dim-HGGvsLGG-nf.png",dpi=200,height = 3,width = 6)

count_dat <- lineage.data %>% group_by(class,lineage.class) %>% summarise(Count=n()) %>% mutate(Frac = Count / sum(Count) * 100) %>% as.data.frame()
count_dat$class <- factor(count_dat$class,levels = c("LGG","HGG"))

ggplot(data=count_dat, mapping=aes(x=lineage.class,y=Frac,fill=class))+
  geom_bar(stat="identity",width=0.5,position='dodge')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='identity',aes(y=Frac+2,label=paste0(round(Frac,0),"%")),color="black",position = position_dodge(0.5),fontface="bold")+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
ggsave("prop_HGG_LGG_2dim_cyc-nf.pdf",width = 5,height = 3)

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=YAP1))+
  geom_point(size=0.1)+
  theme_few()+xlab("AC-like <--------------> NPC-like")+
  scale_color_viridis_c(option = "F")+
  scale_alpha_continuous(range = c(0,1),guide=guide_none())+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(class))
ggsave("2dim-HGGvsLGG-YAP1-nf.png",dpi=200,height = 3,width = 6)
