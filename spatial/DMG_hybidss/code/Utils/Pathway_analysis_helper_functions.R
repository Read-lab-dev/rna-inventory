## Pathway analysis 

## Function to do GO enrichment analysis 
## @para: df (DGE results), 
## @para: df_bg (gives background gene list), 
## @para: dataset (biomart dataset) 
## @para: target cluster ## number 
## @returns: a list of go obj (ego) and go result summary (cluster_summary)
go_analysis <- function(df, df_bg, dataset, cluster){
  ensembl = useMart("ensembl")
  mart <- useDataset(dataset, useMart("ensembl"))
  
  ## Convert ID for df
  merged_gene_ids <- merge(x=df, y=gene_symbol_to_ensembl_id(df$gene, dataset), by.x="gene", by.y="external_gene_name") 
  
  ## Subset significant genes for the target cluster 
  sigOE <- subset(merged_gene_ids, p_val_adj < 0.05)
  sigOE <- sigOE[sigOE$cluster == cluster, ]
  sigOE_genes <- as.character(sigOE$ensembl_gene_id)
  
  ## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
  allOE_genes <- as.character(gene_symbol_to_ensembl_id(rownames(df_bg), dataset)$ensembl_gene_id)
  
  ## Run GO enrichment analysis 
  ego <- enrichGO(gene = sigOE_genes, 
                  universe = allOE_genes, 
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)
  
  print(head(ego))
  cluster_summary <- subset(data.frame(ego), select = -c(geneID))
  return(list(ego=ego, cluster_summary=cluster_summary))
}

go_analysis_easy <- function(genes, bg_genes, dataset = "hsapiens_gene_ensembl", OrgDb = "org.Hs.eg.db"){
  ensembl = useMart("ensembl")
  mart <- useDataset(dataset, useMart("ensembl"))
  
  ## DGE list 
  sigOE_genes <- as.character(gene_symbol_to_ensembl_id(genes, dataset)$ensembl_gene_id)
  
  ## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
  allOE_genes <- as.character(gene_symbol_to_ensembl_id(bg_genes, dataset)$ensembl_gene_id)
  
  ## Run GO enrichment analysis 
  ego <- enrichGO(gene = sigOE_genes, 
                  universe = allOE_genes, 
                  keyType = "ENSEMBL", 
                  OrgDb = OrgDb, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)
  
  print(head(ego))
  cluster_summary <- subset(data.frame(ego), select = -c(geneID))
  return(list(ego=ego, cluster_summary=cluster_summary))
}