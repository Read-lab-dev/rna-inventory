##Subcluster##
setwd("~/rna/BT109VP")
rm(list=ls())
gc()
library(SoupX)
library(Seurat)
library(dplyr)
int.seu <- qs::qread("int.seu.NEW1.qs")
path <- grep(pattern="21047FL-89-01", dir("/home/hzg/rna/seq/ProcessedBam/BT109_VP/"),value = T)
sample_names <- c("Org", "Eng_DMSO", "Eng_VP_S", "Eng_VP_L", "BT109_DMSO", "BT109_VP")

out1 <- load10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-01/outs")
out2 <- load10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-02/outs")
out3 <- load10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-03/outs")
out4 <- load10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-04/outs")
out1 = autoEstCont(out1)
out2 = autoEstCont(out2)
out3 = autoEstCont(out3)
out4 = autoEstCont(out4)
out1 = adjustCounts(out1)
out2 = adjustCounts(out2)
out3 = adjustCounts(out3)
out4 = adjustCounts(out4)

scRNAlist <- list()

# out1 <- Read10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-01/outs/filtered_feature_bc_matrix/")
# out2 <- Read10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-02/outs/filtered_feature_bc_matrix")
# out3 <- Read10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-03/outs/filtered_feature_bc_matrix")
# out4 <- Read10X("/home/hzg/rna/seq/ProcessedBam/BT109_VP/21047FL-89-01-04/outs/filtered_feature_bc_matrix")

for (i in 1:4) {
  scRNAlist[[i]] <- CreateSeuratObject(counts = get(paste0("out",i)),
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
int.seu <- JoinLayers(int.seu)
int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)
int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.8&percent.mt<25&nCount_RNA>1000& 
                    nFeature_RNA >=250)

###SCT##
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
int.seu <- CellCycleScoring(int.seu,s.features = s.genes,g2m.features = g2m.genes)
SCT10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  int.seu <- SCTransform(int.seu,vars.to.regress =  c("percent.mt"),verbose = F)
  int.seu = RunPCA(int.seu,verbose=verbose)
  int.seu <- RunHarmony(int.seu, group.by.vars = "group", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),reduction = "harmony",verbose=verbose)
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),reduction = "harmony",verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
int.seu <- SCT10X(int.seu,res = 0.3)
DimPlot(int.seu,group.by = "celltype",label = T,label.box = T)
FeaturePlot(int.seu,features = "eGFP")

library(harmony)
int.seu$group <- ifelse(int.seu$orig.ident=="Eng_VP_L","Y","N")
int.seu <- NormalizeData(int.seu)
int.seu <- FindVariableFeatures(int.seu,nfeatures = 3000)
int.seu <- ScaleData(int.seu)
int.seu = RunPCA(int.seu,verbose=F)
int.seu <- RunHarmony(int.seu, group.by.vars = "group", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu <- RunTSNE(int.seu,reduction = "harmony",dims = 1:30)
int.seu = FindNeighbors(int.seu,dims=seq(50),verbose=F,reduction = "harmony")
int.seu = FindClusters(int.seu,res=0.4,verbose=F)
FeaturePlot(int.seu,features = "GFAP",order = T)

DimPlot(int.seu,reduction = "umap",label = T)

## DoubletFinder
library(DoubletFinder)
library(Seurat)
library(dplyr)

sweep.res.list <- paramSweep_v3(int.seu, PCs = 1:20, sct = T)
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

library(reticulate)
Sys.setenv(reticulate_python='/home/hzg/miniconda3/envs/r-reticulate/bin/python')

int.seu <- FindClusters(int.seu,resolution = 0.8,algorithm = 4,method = "igraph")

all.markers <- FindAllMarkers(int.seu, only.pos = F, min.pct = 0.25)

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
write.csv(cellmarker,file = "cellmarker_aa.csv")
# FeaturePlot(int.seu,features = c("APOE","GFAP","PDGFRA","MKI67","TP53"),order = T)
DimPlot(int.seu,label = T,label.box = T)

library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",8),8)
DimPlot(int.seu, group.by = "celltype",split.by = "orig.ident",reduction = "umap", label = T,repel = T,label.size = 4,label.box = F,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)&NoAxes()
ggsave("umap_new.png",width = 9,height = 7)

FeaturePlot(int.seu,features = "EPHA3",reduction = "tsne",order = T,split.by = "orig.ident")+theme(legend.position = "right")&NoAxes()


int.seu <- qs::qread("int.seu.NEW.qs")
"#924099" "#DF9ED4" "#406E89" "#544799" "#D2C564" "#D77B5A" "#5DA373" "#CF5A79"
aa <- FindAllMarkers(int.seu,logfc.threshold = 2,only.pos = T)

write.csv(aa,file = "cellmarkers.new.csv")
qs::qsave(int.seu,file = "int.seu.NEW1.qs")

top10.markers <- all.markers[!grepl(c("^ENSG"),rownames(all.markers)),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = rev(p_val_adj))

int.seu$orig.ident <- factor(int.seu$orig.ident,levels = c("Org","Eng_DMSO","Eng_VP_S","Eng_VP_L"))

cell.prop <- as.data.frame(prop.table(table(int.seu$celltype,int.seu$orig.ident),margin = 2))
colnames(cell.prop)<- c("celltype","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=celltype))+
  scale_fill_manual(values = seurat_colors)+
  geom_bar(stat="identity",colour="#222222",width = 0.7)+
  theme_classic()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        panel.border = element_rect(fill = NA,color = "black",linetype = "solid"),
        title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
ggsave("Org_PROP_Celltpye.pdf",dpi = 300,width = 4,height = 4)

DimPlot(int.seu)

