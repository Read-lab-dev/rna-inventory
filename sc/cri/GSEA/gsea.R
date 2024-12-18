.libPaths(c("/home/hzg/R/x86_64-pc-linux-gnu-library/SeuratV5/","/home/hzg/R/x86_64-pc-linux-gnu-library/4.4/"))
library(dplyr)
library(clusterProfiler)

Idents(int.seu) <- int.seu$nmf_cluster
ranked_gene_list <- FindAllMarkers(int.seu)
ranked_gene_lists <- list()
for (i in seq_along(gene_matrix)) {
 ranked_gene <- filter(ranked_gene_list,cluster==i) %>% filter(p_val_adj<0.05) %>% arrange(desc(avg_log2FC))
 ranked_genes <- ranked_gene$avg_log2FC
 names(ranked_genes) <- ranked_gene$gene
 ranked_gene_lists[[colnames(gene_matrix)[i]]] <- ranked_genes
}
ranked_gene_lists <- list()
# Or
# Define the metric for ranking (e.g., log fold change)
# Here, we'll assume the matrix is already log-transformed
# Loop over each column (sample/condition) in the gene matrix
gene_matrix <- as.data.frame(res.rank5@fit@W)
for (i in seq(ncol(gene_matrix))) {
  # Extract the column as a vector
  gene_vector <- gene_matrix[, i]
  names(gene_vector) <- rownames(gene_matrix)
  # Rank the genes based on their expression values
  # Use `order` to sort in decreasing order
  ranked_genes <- sort(gene_vector, decreasing = TRUE)
  # Store the ranked list in the list with the sample name as key
  ranked_gene_lists[[colnames(gene_matrix)[i]]] <- ranked_genes
}

combined_gsea_result <- list()
for (i in colnames(gene_matrix)){
  gsea_input <- ranked_gene_lists[[i]]
  ##Select your Gene Set
  dir='/home/hzg/rna/Bulk_Analysis/MsigDB/'
  gmts <- list.files(dir,pattern = 'gmt')
  gmts
  #Start GSEA Analysis
  gsea_results <- lapply(gmts, function(gmtfile){
    # gmtfile=gmts[2]
    geneset <- read.gmt(file.path(dir,gmtfile)) 
    print(paste0('Now process the ',gmtfile))
    egmt <- GSEA(gsea_input, TERM2GENE=geneset, verbose=FALSE)
    head(egmt)
    return(egmt)
  })
  # gsea_results[[1]] <- setReadable(gsea_results[[1]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  # gsea_results[[2]] <- setReadable(gsea_results[[2]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
  combined_gsea_result[[i]] <- gsea_results_df
}
