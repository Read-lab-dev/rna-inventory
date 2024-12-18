# library(devtools)
# devtools::install_github("yu-lab-vt/CMC@CMC-GOTermActivity")
detach(package:Gotana)
library(Gotana)
library(Seurat)
library(matrixStats)
library(dplyr)
setwd("~/rna/BT109VP")
# Read Raw count
dat <- as.data.frame(GetAssayData(int.seu,slot = "counts"))

dat <- dat %>% dplyr::mutate(Gene=rownames(dat),.before=1)

data.table::fwrite(dat,"Go_data.txt",sep = "\t")


GO_analysis_example <- Read_scData("Go_data.txt", Str_mt = "MT-")

# Quality Control
GO_analysis_example <- QC_scData(GO_analysis_example, Gene_threshold = 0.05,
                                 Count_threshold = 3.5, Cell_threshold = 1500, MT_threshold = 0.2)# Map to GO Term Dataset

GO_analysis_example$AfterQC_List$Gene <- stringr::str_to_title(GO_analysis_example$AfterQC_List$Gene)

GO_analysis_example <- Map_GOSet(GO_analysis_example)

# Run CMC Model
GO_analysis_example <- Run_CMC(GO_analysis_example)

# Calculate GO Term Activity Scores
GO_analysis_example <- GO_Scores(GO_analysis_example)

# Feature Selection: Using Order Statistics Tests
GO_analysis_example <- GO_selected_Order_Statistics(GO_analysis_example)

GO_analysis_example <- GO_analysis_example
# Run PCA
GO_analysis_example <- Run_PCA(GO_analysis_example, Selected_GOs = -1, Score_type = "p_value")

# Build Seurat Object
GO_analysis_example <- creatseurat(GO_analysis_example)

#Analysis using "Seurat" function
Seurat_GO_example <- GO_analysis_example$Seurat_GO

DimHeatmap(Seurat_GO_example, dims = 1, cells = 1000, balanced = TRUE)

Seurat_GO_example <- FindNeighbors(Seurat_GO_example, dims = 1:10)

Seurat_GO_example <- FindClusters(Seurat_GO_example, resolution = 0.5)

Seurat_GO_example <- RunUMAP(Seurat_GO_example, dims = 1:10)

Seurat_GO_example$treat <- substr(colnames(Seurat_GO_example),7,8)
DimPlot(Seurat_GO_example, reduction = "umap", pt.size = 1)
DimPlot(Seurat_GO_example, reduction = "umap", pt.size = 1,group.by = c("seurat_clusters","treat"),label = T)

Idents(Seurat_GO_example) <- Seurat_GO_example$seurat_clusters
cluster.markers <- FindAllMarkers(Seurat_GO_example,min.pct = 0.25)

cluster.markers$ID <- cluster.markers$gene
# cluster.markers$ID <- rownames(cluster.markers)
cluster.markers <- left_join(cluster.markers,GO_analysis_example$GO_Dataset$GO_Term_list,multiple = "any")
top10.markers <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

save(cluster.markers,top10.markers,Seurat_GO_example,file = "Gotana.nsp.Rdata")

# For Further Analysis
GO_analysis_example <- Update_Seurate(GO_analysis_example, Seurat_GO_example)

# Detect DEGs
Cluster_Results <- GO_analysis_example$AfterMapping_List$Cluster
Cluster_0 <- GO_analysis_example$AfterQC_List$Cell_ID[Cluster_Results == '0']
Cluster_1 <- GO_analysis_example$AfterQC_List$Cell_ID[Cluster_Results == '1']
Cluster_5 <- GO_analysis_example$AfterQC_List$Cell_ID[Cluster_Results == '5']
Markers_Cluster0 <- Markers_Detection(GO_analysis_example, cell_group1=Cluster_0, cell_group2=Cluster_1, 
                                      Flag_testAllgene = FALSE, min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE)

# Plot Gene Expression Figure
GenePlot(GO_analysis_example, Gene_list = c("C1ql1", "Tmsb4x", "Ppp1r14b", "Rpl13"), Pt_size = 0.75)

GO_Dataset$Gene_list$Gene <- toupper(GO_Dataset$Gene_list$Gene)

save(GO_Dataset,file = "GO_Dataset.Rdata")

creatseurat <- function (List) 
{
  library(Seurat)
  switch(List$PCA$Score_Type, p_value = Data_Matrix <- List$GO_Term_Activity_Scores$p_Value, 
         Z_score = Data_Matrix <- List$GO_Term_Activity_Scores$Z_Score, 
         `-log10_p` = Data_Matrix <- -log10(List$GO_Term_Activity_Scores$p_Value))
  N_TopGO <- length(List$PCA$Variance)
  Seurat_GO <- CreateSeuratObject(counts = Data_Matrix)
  Seurat_GO <- NormalizeData(Seurat_GO, normalization.method = "LogNormalize", 
                             scale.factor = 1e+05)
  if (List$PCA$Score_Type == "Z_score") {
    tmp_matrix <- Data_Matrix
    tmp_matrix[tmp_matrix == -Inf] <- 0
    Min_scores <- min(tmp_matrix)
    Data_Matrix[Data_Matrix == -Inf] <- 5 * Min_scores
    Seurat_GO@assays[["RNA"]]$data <- Data_Matrix
  }
  if (List$PCA$Score_Type == "p_value") {
    Seurat_GO@assays[["RNA"]]$data <- -log10(Data_Matrix)
  }
  if (List$PCA$Score_Type == "log10_p") {
    Seurat_GO@assays[["RNA"]]$data <- Data_Matrix
  }
  Seurat_GO <- FindVariableFeatures(Seurat_GO, selection.method = "vst", 
                                    nfeatures = N_TopGO)
  Seurat_GO <- ScaleData(Seurat_GO)
  Seurat_GO <- RunPCA(Seurat_GO, features = VariableFeatures(object = Seurat_GO))
  Seurat_GO@reductions[["pca"]]@cell.embeddings <- List$PCA$PCA_matrix
  List[["Seurat_GO"]] <- Seurat_GO
  Seurat_Gene <- CreateSeuratObject(counts = List$AfterQC_List$Data)
  Seurat_Gene <- NormalizeData(Seurat_Gene, normalization.method = "LogNormalize", 
                               scale.factor = 1e+05)
  List[["Seurat_Gene"]] <- Seurat_Gene
  List$AfterQC_List[["Data_Normalized"]] <- Seurat_Gene@assays$RNA$data
  return(List)
}
