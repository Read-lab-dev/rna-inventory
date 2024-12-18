##Subcluster##
setwd("~/rna/BT109VP")
rm(list=ls())
gc()
library(SoupX)
library(Seurat)
library(dplyr)
path <- grep(pattern="21047FL-89-01", dir("/home/hzg/rna/ProcessedBam/BT109_VP"),value = T)
sample_names <- c("Org", "Eng_DMSO", "Eng_VP_S", "Eng_VP_L", "BT109_DMSO", "BT109_VP")
# Load and QC
sc1 = autoEstCont(load10X("/home/hzg/rna/ProcessedBam/BT109_VP/21047FL-89-01-01/outs"))
sc2 = autoEstCont(load10X("/home/hzg/rna/ProcessedBam/BT109_VP/21047FL-89-01-02/outs"))
sc3 = autoEstCont(load10X("/home/hzg/rna/ProcessedBam/BT109_VP/21047FL-89-01-03/outs"))
sc4 = autoEstCont(load10X("/home/hzg/rna/ProcessedBam/BT109_VP/21047FL-89-01-04/outs"))
sc5 = autoEstCont(load10X("/home/hzg/rna/ProcessedBam/BT109_VP/21047FL-89-01-05/outs"))
sc6 = autoEstCont(load10X("/home/hzg/rna/ProcessedBam/BT109_VP/21047FL-89-01-06/outs"))



out1 = adjustCounts(sc1)
out2 = adjustCounts(sc2)
out3 = adjustCounts(sc3)
out4 = adjustCounts(sc4)
out5 = adjustCounts(sc5)
out6 = adjustCounts(sc6)
scRNAlist <- list()
for (i in 1:6) {
  scRNAlist[[i]] <- CreateSeuratObject(counts = get(paste0("out",i)),
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)

int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)

int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.8&percent.mt<25& 
                    nFeature_RNA >=250&nCount_RNA>=1000)
qs::qsave(int.seu,file = "soupx.qs")

standard10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  int.seu = NormalizeData(dat,verbose=verbose)
  int.seu = FindVariableFeatures(int.seu,selection.method = "vst",nfeatures= 3000,verbose=verbose)
  int.seu = ScaleData(int.seu,verbose=verbose)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}

int.seu <- standard10X(int.seu, nPCs=30, res=0.8)
Idents(int.seu) <- int.seu$seurat_clusters
DimPlot(int.seu)
## DoubletFinder
library(DoubletFinder)
library(Seurat)
library(dplyr)

sweep.res.list <- paramSweep_v3(int.seu, PCs = 1:20, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = 0.039                     # 5000细胞对应的doublets rate是3.9%
homotypic.prop <- modelHomotypic(int.seu$seurat_clusters)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(int.seu)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
int.seu <- doubletFinder_v3(int.seu, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, 
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)


DimPlot(int.seu, reduction = "umap", group.by = "DF.classifications_0.25_0.2_2009")

int.seu<- subset(int.seu,subset = DF.classifications_0.25_0.2_2009 =="Singlet")
int.seu <- standard10X(int.seu, nPCs=20, res=0.8)

all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,max.cells.per.ident = 500,
                              logfc.threshold = 0.25) 
top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  
  cluster <- paste(tmp$gene,collapse = ",")
  
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker24.csv")

Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
celltype_col <- sample(MetBrewer::met.brewer("Klimt",17),17)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = celltype_col,label.color = "grey100",
        pt.size = 0.1)+theme(legend.position = "right")&NoAxes()

cell.prop<-as.data.frame(prop.table(table(aa$egfp_true, aa$seurat_clusters)))
colnames(cell.prop)<- c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  scale_fill_manual(values = celltype_col)+
  geom_bar(stat="identity",position="fill")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
qs::qsave(int.seu,file = "int.qs")

################
library(ggplot2)
celltype_col=MetBrewer::met.brewer("Klimt",13)
int.seu <-qs::qread(file = "intnmf.qs")
p1=DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 3,label.box = T,cols = celltype_col,label.color = "grey100",
        pt.size = 0.1)+theme(legend.position = "right")&NoAxes()&NoLegend()
p2=DimPlot(int.seu,group.by = "orig.ident",reduction = "umap",cols = MetBrewer::met.brewer("Klimt",6), label = T,repel = T,label.size = 3,label.box = T,label.color = "grey100",
        pt.size = 0.1)+theme(legend.position = "right")&NoAxes()&NoLegend()
ggsave(p1+p2,filename="umap1.png",height = 4,width = 9)
DimPlot(int.seu,reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = rep(celltype_col,2),label.color = "grey100",
        pt.size = 0.1)+theme(legend.position = "right")&NoAxes()

colors_scale <- rev(met.brewer("Johnson",100))
p3 <- FeaturePlot(int.seu,features = c("AQP4","GFAP","HOPX","SLC1A3","eGFP","EOMES","TBR1","SATB2","RBFOX3"),order = T,pt.size = 0.1)&NoAxes()&NoLegend()

ggsave(p3,filename="umap_feature1.png",height = 8,width = 9)

mark <- FindAllMarkers(gsc,max.cells.per.ident = 250)
top10.markers <- mark[-c(grep(c("^RP[SL]"),rownames(mark)),
                                grep(c("^MT-"),rownames(mark))),] %>% 
  group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

int.seu <- qs::qread("int.qs")
