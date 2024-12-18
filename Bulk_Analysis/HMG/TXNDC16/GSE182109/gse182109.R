###GSE141383
library(data.table)
library(Seurat)
library(tibble)
library(harmony)
gbm.data <- Read10X("./GSE182109/")
int.seu  <- CreateSeuratObject(counts = gbm.data,project = "GSE182109",min.cells=3,min.features=200)

int.seu[["percent.mt"]] <- PercentageFeatureSet(int.seu, pattern = "^MT-")
VlnPlot(int.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
int.seu <- subset(int.seu, subset = percent.mt<25& 
                    nFeature_RNA >=250&nCount_RNA>=1000)
int.seu <- NormalizeData(int.seu)
int.seu <- FindVariableFeatures(int.seu, selection.method = "vst", nfeatures = 3000)
int.seu <- ScaleData(int.seu)
int.seu <- RunPCA(int.seu, verbose = FALSE)
ElbowPlot(int.seu)
int.seu <- RunHarmony(int.seu, group.by.vars = "orig.ident", max_iter = 20,
                      reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
int.seu <- RunUMAP(int.seu,reduction = "harmony",dims = 1:30)
int.seu = FindNeighbors(int.seu,dims=seq(30),verbose=F,reduction = "harmony")
int.seu = FindClusters(int.seu,res=0.2,verbose=F)
DimPlot(int.seu,raster = F,label = T)
FeaturePlot(int.seu,features = "TMEM119",raster = F,order = T)

markers <-FindAllMarkers(int.seu,logfc.threshold = 1,max.cells.per.ident=1000)
top10.markers <- markers[-c(grep(c("^RP[SL]"),rownames(markers)),
                            grep(c("^MT-"),rownames(markers))),] %>% 
  dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
DimPlot(int.seu,label = T)
genes_to_check = c("PTPRC","EPCAM",'PECAM1','MME',"CD3G","CD3E","CD4","CD8A", "CD79A")
DotPlot(int.seu, features = genes_to_check,assay='RNA')  
cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  
  cluster <- paste(tmp$gene,collapse = ",")
  
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker.csv")
qs::qsave(int.seu,file = "./GSE182109.seu.qs")
library(ggplot2)
Idents(int.seu) <- int.seu$seurat_clusters
new.cluster.ids <- clipr::read_clip()##Copy from your handmade-Celltype
names(new.cluster.ids) <- levels(int.seu)
int.seu <- RenameIdents(int.seu, new.cluster.ids)
int.seu$celltype <- Idents(int.seu)
new.levels <- sort(unique(new.cluster.ids))
int.seu$celltype <- factor(int.seu$celltype,levels = new.levels)
Idents(int.seu) <- int.seu$celltype
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",20),20)
DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = seurat_colors,label.color = "grey100",
        pt.size = 0.1,raster = F)
ggsave("UMAP.png",dpi = 200,height = 6)
qs::qsave(int.seu,file = "GSE182109.qs")


mp.seu <- subset(int.seu,celltype=="MP/MG")
mp.seu <- NormalizeData(mp.seu) %>% 
          FindVariableFeatures(nfeatures=3000) %>% 
          ScaleData() %>%
          RunPCA()
mp.seu <- RunHarmony(mp.seu, group.by.vars = "orig.ident", max_iter = 20,
                      reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
ElbowPlot(mp.seu)
mp.seu <- FindNeighbors(mp.seu,reduction = "harmony",dims = 1:30) %>% 
          RunUMAP(reduction = "harmony",dims = 1:30) %>% 
          FindClusters(resolution = 0.5)
DimPlot(mp.seu)
FeaturePlot(mp.seu,features = "TXNDC16",label = T)
markers <-FindAllMarkers(mp.seu,logfc.threshold = 1,max.cells.per.ident=1000)
top10.markers <- markers[-c(grep(c("^RP[SL]"),rownames(markers)),
                            grep(c("^MT-"),rownames(markers))),] %>% 
  dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
