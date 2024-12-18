library(singlet)
library(Seurat)
library(dplyr)
library(cowplot)
set.seed(123)
seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
int.seu <-qs::qread(file = "int.qs")

gsc.seu <- subset(int.seu,idents= c("cGSC","GSC_AS","GSC_OPC"))
gsc.seu <- standard10X(gsc.seu,nPCs = 20)

gsc.seu <- qs::qread("gsc.seu.qs")

all.markers <- FindAllMarkers(gsc.seu, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25) 
top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

gsc.seu <- Seurat::CellCycleScoring(gsc.seu,s.features = cc.genes.updated.2019$s.genes,g2m.features = cc.genes.updated.2019$g2m.genes)
DimPlot(gsc.seu,group.by = "orig.ident",label = T,cols = seurat_colors)+DimPlot(gsc.seu,label = T,cols = rep(seurat_colors,2))

library(NMF)
aa <- NormalizeData(gsc.seu) %>% FindVariableFeatures() %>% ScaleData(do.center = F)
Idents(aa) <- aa$orig.ident 
vm <- gsc.seu@assays$RNA@scale.data
s.f=1:5
res.rank <- nmf(vm, 
                rank = 5,
                seed = "random",
                nrun=5,
                .opt="vp10",
                method = "lee")
qs::qsave(res.rank,file = "res.rank.qs",nthreads = 10L)
plot(res.rank)
fs <- extractFeatures(res.rank, 50L)
fs <- lapply(fs, function(x) rownames(res.rank)[x])
fs <- do.call("rbind",fs)
rownames(fs) <- paste0("cluster", 1:5)
DT::datatable(t(fs))
gsc.seu$nmf_cluster<- predict(res.rank)
Idents(gsc.seu) <- gsc.seu$nmf_cluster
new.cluster.ids <- read.table('cell.txt', sep="\t", header=F)$V1  ##Copy from your handmade-Celltype
Idents(gsc.seu) <- gsc.seu$nmf_cluster
names(new.cluster.ids) <- 1:5
gsc.seu <- RenameIdents(gsc.seu, new.cluster.ids)
gsc.seu$cluster <- Idents(gsc.seu)
gsc.seu$cluster  <- factor(gsc.seu$cluster,levels=new.cluster.ids[c(1,5,2,3,4)])
Idents(gsc.seu)<- gsc.seu$cluster

library(ggplot2)
seurat_prop <- function(obj=int.seu,x="cluster",y="orig.ident",cols=seurat_colors){
cell.prop<-as.data.frame(prop.table(table(obj@meta.data[,x], obj@meta.data[,y])))
colnames(cell.prop)<- c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  scale_fill_manual(values = cols)+
  geom_bar(stat="identity",position="fill")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
}
seurat_prop(obj=gsc.old)+DimPlot(gsc.old,split.by = "orig.ident",cols = seurat_colors)

qs::qsave(gsc.seu,file = "gsc.seu.qs")
gsc.seu$orig_cluster <- paste0(gsc.seu$cluster,"_",gsc.seu$orig.ident)
Idents(gsc.seu)<- gsc.seu$orig_cluster
markers <- FindMarkers(gsc.seu,ident.1 = "eng_VP1",ident.2 = "eng_DMSO" ,test.use = "DESeq2")
int.seu <- qs::qread("int.qs")

apply(NMF::coefficients(res.rank)[s.f,], 2, max)
