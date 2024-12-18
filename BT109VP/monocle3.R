rm(list = ls())
library(monocle3)
library(Seurat)
setwd("~/rna/BT109VP/")
int.seu <- qs::qread("gsc.seu.qs")
int.seu <- subset(int.seu,orig.ident=="Org")
int.seu <- subset(int.seu,celltype=="GSC",invert=T)
expression_data <- as(as.matrix(int.seu@assays$RNA@counts),"sparseMatrix")
cell_data <- int.seu@meta.data
gene_data <- data.frame(gene_short_name = row.names(int.seu), row.names = row.names(int.seu))
cds <- new_cell_data_set(expression_data = expression_data,
                         cell_metadata = cell_data,
                         gene_metadata = gene_data)
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 5)
## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds,cores = 10)
plot_pc_variance_explained(cds)
##Import UMAP axis
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- int.seu@reductions$dim2@cell.embeddings
colnames(int.embed)<- colnames(cds.embed)
cds@int_colData$reducedDims$UMAP <- int.embed
## Step 4: Cluster the cells
cds <- cluster_cells(cds,cluster_method = "leiden")
## Step 5: Learn a graph
cds <- learn_graph(cds,use_partition=F,learn_graph_control=list(minimal_branch_len=5))
plot_cells(cds, color_cells_by = "partition")
plot_cells(cds, color_cells_by = "cluster")
## Step 6: Order cells
cds <- order_cells(cds)#

cds <- reduce_dimension(cds,reduction_method = "LSI",preprocess_method = "LSI")
plot_cells(cds,
           cell_size = 0.8,
           color_cells_by = "pseudotime",
           label_leaves=T,
           label_branch_points=T)+scale_color_viridis_c(option = "I")
ggplot2::ggsave("./monocle/pseudotime_2dim.pdf",height = 5,width = 6)
#
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

#
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

Track_genes_sig <- c("DLG2","TOP2A","SLC1A3","SOX2","NEUROD1","VIM","OPCML","PTPRD","ERBB4")
#
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="celltype", 
                         min_expr=0.5, ncol = 2)
ggplot2::ggsave("./monocle/Genes_Featureplot.png", width = 12, height = 8)
#FeaturePlot
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
ggplot2::ggsave("./monocle/Genes_Featureplot2.png", width = 12, height = 10)

