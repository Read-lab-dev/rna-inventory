# Load libraries

## Libraries for QC and clustering
library(ggplot2)
library(reshape2)
library(dplyr)
library(biomaRt)
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(pagoda2)
library(conos)

## Libraries for additional plots for 
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)


# Define a set of functions 

multiple_QC_metric_plot <- function(qc, metric, metric_short_name, type = c("violin", "histogram", "dot", "bar"), fig_title, ylabel){
  ## Make subsets of df with target variables and melt 
  qc_df = data.frame(qc[,metric])
  colnames(qc_df) = metric_short_name
  qc_metric = melt(qc_df)
  
  if (type == "violin"){
    g = ggplot(qc_metric, aes_string(x="variable", y="value", fill ="variable")) + 
      geom_violin() + ggtitle(fig_title) + ylab(ylabel) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  if (type == "histogram"){
    g = ggplot(qc_metric, aes_string(x="value", fill = "variable")) + 
      geom_histogram() + ggtitle(fig_title) + ylab(ylabel)  +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  if (type == "dot"){
    g = ggplot(qc_metric, aes_string(x="variable", y="value")) + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + ylab(ylabel) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  if (type == "bar"){
    g = ggplot(qc_metric, aes_string(x="variable", y="value", color="variable")) + stat_summary(fun.data="mean_sdl", mult=0.5, geom="crossbar", width=0.1) + ylab(ylabel) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  g
}

## Function to make QC comparison between two samples 
## @para: qc1 and qc2, two qc datasets
## @para: metric, a vector of qc metrics to compare 
## @para: metric_short_name: a vector of qc metrics to be labeled on plots (one to one correspond to metric)
## @para: comparison, samples to be compared (default to "Cell" vs "Nucleus")
## @para: type, "violin" or "histogram" or "dot" or "bar" graphs
## @para: fig_title, title of plots
## @para: ylabel, label of y axis in the plots
## @return: a plot that compares a QC metric 
QC_metric_plots_dual_samples <- function(qc1, qc2, metric, metric_short_name, comparison = c("Cell", "Nucleus"), type = c("violin", "histogram", "dot", "bar"), fig_title, ylabel){
  ## Make subsets of df with target variables; add Type to each df 
  ## Merge into a single df; melt 
  qc_1 = data.frame(qc1[,metric])
  colnames(qc_1) = metric_short_name
  qc_1$Type = comparison[1]
  qc_2 = data.frame(qc2[,metric])
  colnames(qc_2) = metric_short_name
  qc_2$Type = comparison[2]
  qc_metric = rbind(qc_1,qc_2)
  qc_metric = melt(qc_metric)
  
  ggplot(qc_metric, aes_string(x="Type", y="value", color="Type")) + 
    geom_violin() + facet_grid(~variable) + ggtitle(fig_title) + ylab(ylabel) + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  
  if (type == "violin"){
    g = ggplot(qc_metric, aes_string(x="Type", y="value", fill ="Type")) + 
      geom_violin() + facet_grid(~variable) + ggtitle(fig_title) + ylab(ylabel) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  if (type == "histogram"){
    g = ggplot(qc_metric, aes_string(x="value", fill = "Type")) + 
      geom_histogram() + facet_grid(~variable) + ggtitle(fig_title) + ylab(ylabel)  +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  if (type == "dot"){
    g = ggplot(qc_metric, aes_string(x="Type", y="value")) + geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + facet_grid(~variable) + ylab(ylabel) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  if (type == "bar"){
    g = ggplot(qc_metric, aes_string(x="Type", y="value", color="Type")) + stat_summary(fun.data="mean_sdl", mult=0.5, geom="crossbar", width=0.1) +  facet_grid(~variable) + ylab(ylabel) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  g
}

## Function to make gene_id as rownames for a cm
## @para: df (a raw cm still with a gene_id column)
## @returns: df with gene_id removed, and rownames as gene_id 
addRownames <- function(df){
  ensembl_id = sapply(df$gene_id, function(x) unlist(strsplit(x, "\\."))[1])
  df$ensembl_id = ensembl_id
  df = distinct(df, ensembl_id, .keep_all = TRUE)
  rownames(df) = df$ensembl_id
  df = subset(df, select=-c(gene_id, ensembl_id))
}

## Function to convert gene_id into gene symbols
## @para: df, a cm still with gene_id as rownames
## @para: dataset, biomart dataset ("hsapiens_gene_ensembl" or "mmusculus_gene_ensembl")
## @returns: df with gene symbol as rownames 
ensembl_to_symbol_biomart <- function(df, dataset = c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl")){
  ensembl = useMart("ensembl")
  mart <- useDataset(dataset, useMart("ensembl"))
  ids <- rownames(df)
  genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name"), values=ids, mart=mart)
  genes <- genes[which(duplicated(genes$external_gene_name) == FALSE), ] 
  merged_gene_ids <- merge(x=df, y=genes, by.x="row.names", by.y="ensembl_gene_id")  
  rownames(merged_gene_ids) = merged_gene_ids$external_gene_name
  merged_gene_ids = subset(merged_gene_ids, select=-c(Row.names, external_gene_name))
}

## Function to convert gene_id into gene symbols
## @para: df, a cm still with gene_id as rownames
## @para: dataset, biomart dataset ("hsapiens_gene_ensembl" or "mmusculus_gene_ensembl")
## @returns: df with gene symbol as rownames 
symbol_to_entrez_biomart <- function(df, dataset = c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl")){
  ensembl = useMart("ensembl")
  mart <- useDataset(dataset, useMart("ensembl"))
  ids <- rownames(df)
  genes <- getBM(filters="external_gene_name", attributes=c("external_gene_name", "entrezgene"), values=ids, mart=mart)
  genes <- genes[which(duplicated(genes$entrezgene) == FALSE), ] 
  merged_gene_ids <- merge(x=df, y=genes, by.x="row.names", by.y="external_gene_name")  
  merged_gene_ids <- na.omit(merged_gene_ids)
  rownames(merged_gene_ids) = merged_gene_ids$entrezgene
  merged_gene_ids = subset(merged_gene_ids, select=-c(Row.names, entrezgene))
}

## Create a scater object 
## @para: count_matrix, cm
## @para: mito_genes, a list of mitocondrial genes as gene symbols 
## @return: a scater object 
createScaterObj <- function(count_matrix, mito_genes){
  ## create obj
  obj = SingleCellExperiment(assays = list(counts = as.matrix(count_matrix)))
  ## Add mitocondrial gene controls 
  isSpike(obj, "MT") <- rownames(count_matrix) %in% mito_genes
  ## Calculate QC metrics
  obj <- calculateQCMetrics(
    obj,
    feature_controls = list(
      MT = isSpike(obj, "MT")
    )
  )
  return(obj)
}

## Construct a seurat object
## @para: count_matrix, cm
## @para: mito_symbol, mitocondrial genes in gene name/symbol format
## @para: min.cells, minimal number of cells that express a certain gene for that gene to be considered as expressed (default = 1)
## @para: is.expr, gene expression detection limit; default = 0 (for UMI); reads should use 5
## @para: project, project name for this seurat object 
seuratObj <- function(count_matrix, mito_symbol, min.cells = 1, is.expr = 0, project){
  seurat_obj <- CreateSeuratObject(round(count_matrix), min.cells = min.cells, is.expr = is.expr, project = project)
  percent.mito <- Matrix::colSums(seurat_obj@raw.data[mito_symbol, ]) / Matrix::colSums(seurat_obj@raw.data)
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.mito, col.name = "percent.mito")
}

multiple_seurat_QC_metric_plot_single_sample <- function(seurat_obj, metric, metric_short_name, log_tran = TRUE, plot_type = c("violin", "histogram"), fig_title, xlabel, ylabel){
  ## Combine the target attributes of each dataset into one and melt it 
  qc = data.frame(seurat_obj@meta.data[,metric])
  colnames(qc) = metric_short_name
  if (log_tran == TRUE){
    qc = log10(qc+1)
  }
  qc_metric = melt(qc)
  
  ## Plot violin or histogram 
  if (plot_type == "violin"){
    g <- ggplot(qc_metric, aes_string(x="variable", y="value", fill="variable")) + 
      geom_violin() + facet_grid(.~variable) + ggtitle(fig_title) + xlab(xlabel) + ylab(ylabel) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  
  if (plot_type == "histogram"){
    g <- ggplot(qc_metric, aes_string("value", fill="variable")) + 
      geom_histogram() + facet_grid(.~variable) + ggtitle(fig_title) + xlab(xlabel) + ylab(ylabel) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  
  g
}

## Plot total number of reads/genes using seurat object
## @para: seurat_obj1/2, a seurat object for each type to be compared 
## @para: metric, a vector of qc metrics to compare 
## @para: metric_short_name: a vector of qc metrics to be labeled on plots (one to one correspond to metric)
## @log_tran: whether log10 transform the data
## @para: plot_type, violin or histogram
## @para: fig_title, title of plots
## @para: ylabel, label of y axis in the plots
## @return: a plot that compares a QC metric 
multiple_seurat_QC_metric_plot <- function(seurat_obj1, seurat_obj2, metric, metric_short_name, comparison = c("Cell", "Nucleus"), log_tran = TRUE, plot_type = c("violin", "histogram"), fig_title, xlabel, ylabel){
  ## Combine the target attributes of each dataset into one and melt it 
  qc1 = data.frame(seurat_obj1@meta.data[,metric])
  colnames(qc1) = metric_short_name
  if (log_tran == TRUE){
    qc1 = log10(qc1+1)
  }
  qc1$Type = comparison[1]
  
  qc2 = data.frame(seurat_obj2@meta.data[,metric])
  colnames(qc2) = metric_short_name
  if (log_tran == T){
    qc2 = log10(qc2+1)
  }
  qc2$Type = comparison[2]
  
  qc_metric = rbind(qc1,qc2)
  qc_metric = melt(qc_metric)
  
  ## Plot violin or histogram 
  if (plot_type == "violin"){
    g <- ggplot(qc_metric, aes_string(x="Type", y="value", fill="Type")) + 
      geom_violin() + facet_grid(.~variable) + ggtitle(fig_title) + xlab(xlabel) + ylab(ylabel) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  
  if (plot_type == "histogram"){
    g <- ggplot(qc_metric, aes_string("value", fill="Type")) + 
      geom_histogram() + facet_grid(Type~variable) + ggtitle(fig_title) + xlab(xlabel) + ylab(ylabel) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
  }
  
  g
}

## Compute % of reads mapped to each gene type 
## @para: count_matrix, cm
## @para: dataset, biomart dataset to use 
## @return: a df with % of each gene type for each cell
pct_gene_type <- function(count_matrix, dataset = c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl")){
  toMatch <- c("*pseudogene", "antisense", "*lincRNA", "miRNA", "misc_RNA", 
               "protein_coding", "snRNA", "snoRNA")
  gene_type = NULL
  ensembl = useMart("ensembl")
  mart <- useDataset(dataset, useMart("ensembl"))
  
  for (i in 1:dim(count_matrix)[2]){
    print(colnames(count_matrix)[i])
    tmp = count_matrix[,i]
    names(tmp) = rownames(count_matrix)
    tmp = tmp[tmp>0]
    gene_biotypes = getBM(attributes=c("external_gene_name", "gene_biotype"), filters='external_gene_name', values=c(names(tmp)), mart=mart)$gene_biotype
    gene_biotypes_dist = NULL
    for (term in toMatch){
      gene_biotypes_dist = c(gene_biotypes_dist, length(grep(term, gene_biotypes)))
    }
    gene_biotypes_dist = gene_biotypes_dist/sum(gene_biotypes_dist) * 100
    gene_type = rbind(gene_type, gene_biotypes_dist)
  }
  
  gene_type = data.frame(gene_type)
  colnames(gene_type) = c("pseudogene", "antisense", "lincRNA", "miRNA", "misc_RNA", 
                          "protein_coding", "snRNA", "snoRNA")
  gene_type$cell_id = colnames(count_matrix)
  gene_type = melt(gene_type)
  colnames(gene_type) = c("cell_id", "gene_type", "pct")
  return(gene_type)
}

## Plot gene type 
## @para: df_gene_type, df with gene type %
## @para: title, title of the plot 
## @para: type, "(Cell)" or "(Nucleus)"
## @para: xlabel, x axis label 
## @para: ylabel, y axis label 
## @returns: a plot 
plot_gene_type <- function(df_gene_type, title="Gene types detected", type = c("(Cell)", "(Nucleus)"), xlabel="Cell_ID", ylabel="Pct"){
  plot_title = paste(title,type, sep = " ")
  ggplot(df_gene_type, aes(x=cell_id, y=pct, fill=gene_type)) + geom_bar(stat="identity") + scale_fill_brewer(palette="Dark2") + ggtitle(plot_title) + xlab(xlabel) + ylab(ylabel) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=20), axis.text.x=element_blank())
}

## Function to convert a list of gene symbols to a df with gene symbol and corresponding ensembl id
## @para: gene (a list of gene symbols)
## @para: dataset (biomart dataset)
## @returns: df with matched ensembl_gene_id and external_gene_name 
gene_symbol_to_ensembl_id <- function(gene, dataset){
  ensembl = useMart("ensembl")
  mart <- useDataset(dataset, useMart("ensembl"))
  ids <- getBM(filters="external_gene_name", attributes=c("ensembl_gene_id", "external_gene_name"), values=gene, mart=mart)
  non_duplicates <- which(duplicated(ids$external_gene_name) == FALSE)
  ids <- ids[non_duplicates, ] 
}

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

preprocessing_ss2_data <- function(cm, total_gene_cutoff = 500, gene_expr = 5, gene_per_cell = 2){
  ## Filter out cells with total number of genes < 500
  cm <- cm[,which(apply(cm, 2, function(x) sum(x > gene_expr) > total_gene_cutoff))]
  ## Filter out genes that are expressed in 0 or 1 cell 
  filter_genes <- apply(cm, 1, function(x) length(x[x > gene_expr]) >= gene_per_cell)
  cm <- cm[filter_genes,]
  ## Filter out duplicated gene names (usually do not exist after converting ensembl-ID to gene symbols)
  rownames(cm) <- make.unique(rownames(cm))
  ## Convert df to matrix and return
  cm <- as.matrix(cm)
  return(cm)
}