##Subcluster##
setwd("~/rna/sc/cri")
rm(list=ls())
gc()
library(SoupX)
library(Seurat)
library(dplyr)
path <- dir("./raw/")
sample_names <- c("DMSO", "cz_long", "cz_short")
scRNAlist <- list()
for (i in 1:3) {
  testdf <- Read10X(paste0("./raw/",path[i]))
  scRNAlist[[i]] <- CreateSeuratObject(counts = testdf,
                                       project = sample_names[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_names[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)

int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
int.seu$log10GenesPerUMI <- log10(int.seu$nFeature_RNA)/log10(int.seu$nCount_RNA)
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)

int.seu <- subset(int.seu, subset = log10GenesPerUMI>0.8&percent.mt<25&nCount_RNA>1000&nCount_RNA<50000& 
                    nFeature_RNA >=250)
options(future.globals.maxSize= 891289600000)
standard10X = function(dat,nPCs=30,res=1.0,verbose=FALSE){
  int.seu = NormalizeData(int.seu)
  int.seu = FindVariableFeatures(int.seu,nfeatures = 3000)
  int.seu = ScaleData(int.seu)
  int.seu = RunPCA(int.seu,verbose=verbose)
  # library(harmony)
  # int.seu = RunHarmony(int.seu, group.by.vars = "orig.ident", reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
  int.seu = RunUMAP(int.seu,dims=seq(nPCs),reduction = "pca",verbose=verbose)
  int.seu = FindNeighbors(int.seu,dims=seq(nPCs),reduction = "pca",verbose=verbose)
  int.seu = FindClusters(int.seu,res=res,verbose=verbose)
  return(int.seu)
}
int.seu <- standard10X(int.seu,res = 0.3)

DimPlot(int.seu,label = T,label.box = T)
FeaturePlot(int.seu,features = "eGFP")
int.seu = FindClusters(int.seu,res=0.3,method = "igraph",algorithm = 4)


## DoubletFinder
library(DoubletFinder)
library(Seurat)
library(dplyr)
int.seu[["RNA"]] <- as(object = int.seu[["RNA"]], Class = "Assay")
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

int.seu$double_cell <- int.seu$DF.classifications_0.25_0.001_1223

int.seu$DF.classifications_0.25_0.001_1223 <- NULL
int.seu <-subset(int.seu,double_cell=="Singlet")
int.seu <- JoinLayers(int.seu)
Idents(int.seu) <- int.seu$seurat_clusters
all.markers <- FindAllMarkers(int.seu, only.pos = T, min.pct = 0.25,logfc.threshold = 1)
write.csv(all.markers,file = "cellmarker_all.csv")

top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^ENSG"),rownames(all.markers)),
                                grep(c("^MIR"),rownames(all.markers)),
                                grep(c("^LINC"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  cluster <- paste(tmp$gene,collapse = ",")
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker.csv")
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
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",11),11)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1)
ggsave("umap_new1.png",width = 9,height = 6)

# FeaturePlot(int.seu,features = "EPHA3",reduction = "tsne",order = T,split.by = "orig.ident")+theme(legend.position = "right")&NoAxes()
# 
# int.seu <- RunPCA(int.seu, assay = "SCT", npcs = 50)
# 
# DimPlot(int.seu,label = T)
# 
int.seu <- qs::qread("int.seu.qs")
# 
# aa <- FindAllMarkers(int.seu,test.use = "MAST",logfc.threshold = 1,only.pos = T)
# 
# write.csv(aa,file = "cellmarkers.new.csv")
# qs::qsave(int.seu,file = "int.seu.qs")
# 
# top10.markers <- all.markers[!grepl(c("^ENSG"),rownames(all.markers)),] %>% 
#   group_by(cluster) %>% top_n(n = 10, wt = rev(p_val_adj))
# 
# int.seu$orig.ident <- factor(int.seu$orig.ident,levels = c("Org","Eng_DMSO","Eng_VP_S","Eng_VP_L"))
# cell.prop <- as.data.frame(prop.table(table(int.seu$celltype,int.seu$orig.ident),margin = 2))
# colnames(cell.prop)<- c("celltype","sample","proportion")
# ggplot(cell.prop,aes(sample,proportion,fill=celltype))+
#   scale_fill_manual(values = seurat_colors)+
#   geom_bar(stat="identity",colour="#222222",width = 0.7)+
#   theme_classic()+
#   theme(axis.ticks.length=unit(0.1,'cm'),
#         panel.border = element_rect(fill = NA,color = "black",linetype = "solid"),
#         title = element_blank(),
#         axis.text.x = element_text(angle = 45,hjust = 1))+
#   guides(fill=guide_legend(title=NULL))
# ggsave("Org_PROP_Celltpye1.pdf",dpi = 300,width = 4,height = 4)

int.seu <- qs::qread("tumor.qs")
int.seu$Barcode <- rownames(int.seu@meta.data)
fusion <- read.table("./FUSION/reads_per_barcode_ZM-fusion.txt",header = T)
fusion <- read.table("./FUSION/reads_per_barcode_ex14.txt",header = T)
fusion$Sample <- gsub("CzL","cz_long",fusion$Sample);fusion$Sample <- gsub("CzS","cz_short",fusion$Sample)
fusion$Barcode <- paste0(fusion$Sample,"_",fusion$Barcode,"-1")
metadata <- dplyr::left_join(int.seu@meta.data,fusion)
int.seu$ZM <- ifelse(is.na(metadata$Count),"not-detected","fusion")
int.seu$ex14 <- ifelse(is.na(metadata$Count),"not-detected","ex14-skipping")
table(int.seu$ZM,int.seu$orig.ident,int.seu$double_cell)
Idents(int.seu) <- int.seu$ZM
DimPlot(int.seu,cells.highlight = WhichCells(int.seu,idents="Fusion"),split.by = "double_cell")+FeaturePlot(int.seu,features = "eGFP",split.by = "double_cell")

DotPlot(int.seu,features = c("MET","PTPRZ1","EGFR","eGFP"))+DimPlot(int.seu,label = T)

DimPlot(int.seu,group.by = c("ZM","ex14"),cols = c("darkred","grey"),reduction = "dim2")
int.seu$ZM <- ifelse(int.seu$ZM)

seurat_colors <- sample(MetBrewer::met.brewer("Klimt",8),8)
cell.prop <- as.data.frame(prop.table(table(int.seu$ex14,int.seu$orig.ident),margin = 2))
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

Idents(int.seu) <- int.seu$orig.ident
VlnPlot(int.seu,features = c("MET"),split.by = c("ZM"),pt.size = 0.01)
qs::qsave(int.seu,file = "tumor.qs")

FeaturePlot(int.seu,features = "MET",reduction = "dim2")
DimPlot(int.seu,cells.highlight = WhichCells(int.seu,expression = MET>0),reduction = "dim2")
