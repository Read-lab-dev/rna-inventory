###GSE141383
library(data.table)
library(Seurat)
library(tibble)
library(harmony)
file <- dir(path = "./GSE141383_RAW/",pattern="*.txt.gz")
filename <- gsub(".cou","",substr(file,12,20))
filename <- gsub(".co","",filename)
scRNAlist <- list()
for (i in 1:length(file)) {
  SCcounts <- fread(paste0("./GSE141383_RAW/",file[i]))[,-1]
  SCcounts <- SCcounts[!duplicated(SCcounts$Gene),]
  SCcounts <- column_to_rownames(SCcounts,var = "Gene")
  scRNAlist[[i]] <- CreateSeuratObject(counts = SCcounts,
                                       project = filename[i],min.cells=3,min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=filename[i])
}
int.seu <- Reduce(function(x,y) merge(x,y) , scRNAlist)
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
int.seu = FindClusters(int.seu,res=0.5,verbose=F)
DimPlot(int.seu)+DimPlot(int.seu,group.by = "orig.ident",label = T)

markers <-FindAllMarkers(int.seu,logfc.threshold = 1)
top10.markers <- markers[-c(grep(c("^RP[SL]"),rownames(markers)),
                            grep(c("^MT-"),rownames(markers))),] %>% 
  dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
DimPlot(int.seu,label = T)
genes_to_check = c("PTPRC","EPCAM",'PECAM1','MME',"CD3G","CD3E", "CD79A")
DotPlot(int.seu, features = genes_to_check,assay='RNA')  
cellmarker <- NULL
for (i in levels(int.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  
  cluster <- paste(tmp$gene,collapse = ",")
  
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(int.seu))
write.csv(cellmarker,file = "cellmarker.csv")
qs::qsave(int.seu,file = "./GSE.seu.qs")

mp.seu <- subset(int.seu,CellAssignment=="Macrophage")
mp.seu <- NormalizeData(mp.seu) %>% 
          FindVariableFeatures(nfeatures=3000) %>% 
          ScaleData() %>%
          RunPCA()
ElbowPlot(mp.seu)
mp.seu <- FindNeighbors(mp.seu,dims = 1:15) %>% 
          RunUMAP(dims=1:15) %>% 
          FindClusters(resolution = 0.5)
DimPlot(mp.seu)
FeaturePlot(mp.seu,features = "CD163")
