#scAB
library(scAB)
library(Seurat)
library(preprocessCore)
data("data_survival")

bulk_dataset <- expr_norm

phenotype <- ifelse(metadata$Treat=="sh3",1,0)

sc_dataset <- run_seurat(gbm.seu,verbose = FALSE)

UMAP_celltype <- DimPlot(sc_dataset,label = T,reduction = "dim2")

scAB_data <- create_scAB(sc_dataset,bulk_dataset,phenotype,method = "binary")

K <- select_K(scAB_data)
K

scAB_result <- scAB(Object=scAB_data, K=K)
sc_dataset <- findSubset(sc_dataset, scAB_Object = scAB_result, tred = 2)

UMAP_scAB <- DimPlot(sc_dataset,reduction = "dim2",group.by="scAB_Subset1",cols=c("#80b1d3","red"),pt.size=0.001,order=c("scAB+ cells","Other cells"))
patchwork::wrap_plots(plots = list(UMAP_celltype,UMAP_scAB), ncol = 2)


