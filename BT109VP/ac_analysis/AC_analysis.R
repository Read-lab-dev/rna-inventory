library(Seurat)
library(ggplot2)
setwd("~/rna/BT109VP/ac_analysis")
int.seu <- qs::qread("../int.seu.NEW1.qs")
DimPlot(int.seu,label = T)+
  DotPlot(int.seu,features = c("SLC1A3","GFAP","AQP4","APOE"))

ac.seu$orig.ident <- factor(ac.seu$orig.ident,levels = c("Org","Eng_DMSO","Eng_VP_S","Eng_VP_L"))

ac.seu <- subset(int.seu,celltype%in%c("AC","RG"))

standard10X = function(dat,nPCs=50,res=0.5,verbose=FALSE){
  int.seu = NormalizeData(dat,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = harmony::RunHarmony(int.seu,group.by.vars="group")
  int.seu = RunUMAP(int.seu,reduction = "harmony",seq(nPCs))
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),verbose=verbose,reduction = "harmony")
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
ac.seu <- standard10X(ac.seu,res = 0.5)
num <- length(levels(Idents(ac.seu)))
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",num),num)
DimPlot(ac.seu,group.by = "orig.ident")+DimPlot(ac.seu)
cell.prop<-as.data.frame(prop.table(table(ac.seu$seurat_clusters, ac.seu$orig.ident)))
colnames(cell.prop)<- c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  scale_fill_manual(values = seurat_colors)+
  geom_bar(stat="identity",position="fill")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))+
DimPlot(ac.seu,label = T,label.box = T,cols = seurat_colors)

all.markers <- FindAllMarkers(ac.seu, only.pos = T, min.pct = 0.25)
top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
cellmarker <- NULL
for (i in levels(ac.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(ac.seu))
write.csv(cellmarker,file = "cellmarker_ac.csv")
qs::qsave(ac.seu,file = "ac.seu.qs")
######Velocyto########
aa <- gsub("-1","x",Cells(ac.seu))
aa <- gsub("Org_","21047FL-89-01-01:",aa)
aa <- gsub("Eng_DMSO_","21047FL-89-01-02:",aa)
aa <- gsub("Eng_VP_S_","21047FL-89-01-03:",aa)
aa <- gsub("Eng_VP_L_","21047FL-89-01-04:",aa)
emb <- Embeddings(ac.seu, reduction = "umap")
rownames(emb) <- aa
write.csv(aa, file = "cellID_obs_ac.csv", row.names = F)
write.csv(emb, file = "cell_embeddings_ac.csv")
bb <- data.frame(CellID=aa,x=Idents(ac.seu))
write.csv(bb, file = "clusters_ac.csv",row.names = F)

#####Reactive markers
new.sig.id <- read.csv("signature.csv",na.strings = "")
ac.seu <- AddModuleScore(ac.seu,nbin = 30,ctrl = 100,features = list(na.omit(new.sig.id$Reactive)),name = "Reactive")
ac.seu <- AddModuleScore(ac.seu,nbin = 30,ctrl = 100,features = list(na.omit(new.sig.id$Fetal)),name = "Fetal")
ac.seu <- AddModuleScore(ac.seu,nbin = 30,ctrl = 100,features = list(na.omit(new.sig.id$Mature)),name = "Mature")

VlnPlot(ac.seu,features = c("Reactive1"),alpha = 0.5)
ggsave("ac-signature.png",width = 15,height = 4)
Idents(ac.seu) <- ac.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(ac.seu)
ac.seu <- RenameIdents(ac.seu, new.cluster.ids)
ac.seu$celltype <- Idents(ac.seu)
new.levels <- sort(unique(new.cluster.ids))
ac.seu$celltype <- factor(ac.seu$celltype,levels = new.levels)
Idents(ac.seu) <- ac.seu$celltype

######Prop in subset#########
ac.seu <- subset(int.seu,celltype%in%c("InN"))
ac.seu <- subset(int.seu,DLG4>0&celltype%in%c("InN","imN","InN3"))
ac.seu <- standard10X(ac.seu,nPCs=15,res=0.5)
cell.prop <- as.data.frame(prop.table(table(ac.seu$orig.ident)))
colnames(cell.prop)<- c("sample","proportion")
cell.prop$sample <- factor(cell.prop$sample,levels = c("Org","Eng_DMSO","Eng_VP_S","Eng_VP_L"))
p1=ggplot(cell.prop,aes(sample,proportion,fill=sample),)+
  scale_fill_manual(values = seurat_colors)+
  geom_bar(stat="identity")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
p2=DimPlot(ac.seu,label = T,cols = seurat_colors)
p1/p2

######Prop in DLG4#########
ac.seu <- subset(int.seu,NEFL>0)
cell.prop <- as.data.frame(prop.table(table(ac.seu$celltype)))
colnames(cell.prop)<- c("sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=sample),)+
  scale_fill_manual(values = seurat_colors)+
  geom_bar(stat="identity")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
DimPlot(ac.seu,label = T,cols = seurat_colors)
p1/p2

######Prop in AQP4#########
ac.seu <- subset(int.seu,AQP4>0)
cell.prop <- as.data.frame(prop.table(table(ac.seu$celltype)))
colnames(cell.prop)<- c("sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=sample),)+
  scale_fill_manual(values = seurat_colors)+
  geom_bar(stat="identity")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
ac.seu <- subset(int.seu,AQP4>0&celltype%in%c("AC","AC2"))

cell.prop <- as.data.frame(prop.table(table(ac.seu$orig.ident)))
ac.seu <- standard10X(ac.seu,nPCs=15,res=0.5)
colnames(cell.prop)<- c("sample","proportion")
cell.prop$sample <- factor(cell.prop$sample,levels = c("Org","Eng_DMSO","Eng_VP_S","Eng_VP_L"))
p1=ggplot(cell.prop,aes(sample,proportion,fill=sample),)+
  scale_fill_manual(values = seurat_colors)+
  geom_bar(stat="identity")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
p2=FeaturePlot(int.seu,features = "AQP4",order = T)
p1/p2

meta_list <- c("ALDOC","FABP7","MAOB","TSPO","SOX9","S100B","SLC1A3","SLC1A2","STAT3")
ac.seu <- AddModuleScore(ac.seu,features = meta_list,nbin = 30,name = "AC")
VlnPlot(ac.seu,features = "AC1",group.by = "orig.ident")

qs::qsave(ac.seu,file = "ac.seu.qs")
aa <- FindAllMarkers(ac.seu,logfc.threshold = 1)
top10.markers <- aa%>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

Idents(ac.seu) <- ac.seu$seurat_clusters
FeaturePlot(ac.seu,features = "Reactive1")

