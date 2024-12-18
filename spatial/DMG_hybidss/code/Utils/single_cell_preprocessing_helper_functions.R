# Load libraries

## Libraries for basic preprocessing 
library(reshape2)
library(dplyr)

## Single cell libraries
library(Seurat)
library(pagoda2)
##library(conos)

## Libraries for plotting 
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

## Annotation and pathway analysis library 
##library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

## Others
library(Rtsne)

# Define functions

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## ==========================================================
## Functions for preprocessing count matrix 
## ==========================================================

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

## Function to convert gene_id into gene symbols using biomart
## @para: df, a cm still with gene_id as rownames
## @para: dataset, biomart dataset ("hsapiens_gene_ensembl" or "mmusculus_gene_ensembl")
## @returns: df with gene symbol as rownames 
ensembl_to_symbol_biomart <- function(df, dataset = c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl")){
  ## Initiate ensembl db
  ensembl = useMart("ensembl")
  mart <- useDataset(dataset, useMart("ensembl"))
  ## Look up HUGO from ensembl gene id (with version) and remove duplicates
  ids <- rownames(df)
  genes <- getBM(filters="ensembl_gene_id_version", attributes=c("ensembl_gene_id_version", "external_gene_name"), values=ids, mart=mart)
  genes <- genes[which(duplicated(genes$external_gene_name) == FALSE), ] 
  ##merged_gene_ids <- merge(x=df, y=genes, by.x="row.names", by.y="ensembl_gene_id")  
  ##rownames(merged_gene_ids) = merged_gene_ids$external_gene_name
  ##merged_gene_ids = subset(merged_gene_ids, select=-c(Row.names, external_gene_name))
  ## Make a vector of HUGO, named by ensembl id 
  gene_symbols = genes$external_gene_name
  names(gene_symbols) = genes$ensembl_gene_id_version
  ## Remove na and duplicates
  gene_symbols = gene_symbols[!is.na(gene_symbols)]
  gene_symbols = gene_symbols[!duplicated(gene_symbols)]
  ## Convert rownames of df (cm) from ensembl id to HUGO 
  df = df[names(gene_symbols),]
  rownames(df) = gene_symbols
  return(df)
}

## Function to convert gene_id into gene symbols using AnnotationDB
## @para: cm, a cm still with gene_id as rownames
ensembl_to_symbol_annotationdb <- function(cm){
  ## Remove duplicated gene ids if version number is removed 
  gene_names = sapply(rownames(cm), function(x) unlist(strsplit(x, split = ".", fixed=T))[1])
  names(gene_names) = rownames(cm)
  gene_names = gene_names[!duplicated(gene_names)]
  ## Convert gene_names from gene ids into gene symbols 
  cm = cm[names(gene_names),]
  rownames(cm) = gene_names
  gene_symbols = mapIds(org.Hs.eg.db, keys = rownames(cm), column = "SYMBOL", multiVals = "first", keytype = "ENSEMBL") 
  ## Remove NA and duplicates from gene_symbols and make them rownames of cm
  gene_symbols = gene_symbols[!is.na(gene_symbols)]
  gene_symbols = gene_symbols[!duplicated(gene_symbols)]
  cm = cm[names(gene_symbols),]
  rownames(cm) = gene_symbols
  return(cm)
}

## Remove "_" and "." from cell names
## @para: cell_name, a vector of cell names 
cleanCellName <- function(cell_name){
  cell_name = gsub("_", "", cell_name, fixed = T)
  cell_name = gsub("-", ".", cell_name, fixed = T)
  return(cell_name)
}

## ==========================================================
## Functions for QC and filter smart-seq2 matrix 
## ==========================================================

## normalize and filtering tpm dataset 
## @param cm count matrix
## @param filter_cell whether to filter cells, default = TRUE
## @param scale_factor scaling factor to shrink tpm, default = 10
## @param gene_detect gene detection limit, default = 1
## @param min_genes minimal number of genes for calling a cell, default = 2500
## @param min_mean_hg minimal mean expression of house keeping genes for calling a cell, default = 2.5
## @param hg_list the list of house keeping genes, default = 
## @param gene_cutoff average expression cutoff to keep a gene 
## @param centering whether to center filtered gene expressions, default = TRUE
## @param verbose whether to print summary statistics, default = TRUE
normalizeTPM <- function(cm, filter_cell = TRUE, log_base=2, scale_factor=10, 
                         gene_detect=1, min_genes=2500, 
                          min_mean_hg=2.5, hg_list=NULL, 
                         filter_gene_relaxed=F, gene_cutoff=4, gene_cell_cutoff=10, verbose = TRUE){
  cm = as.matrix(cm)
  message("The number of cells in the unfilered cm:", dim(cm)[2])
  if (filter_cell){
    ## Filter cells with low total number of genes detected
    message("Filtering based on total number of detected genes")
    cm = cm[,which(colSums(cm >= gene_detect) >= min_genes)]
    message("Number of cells left after total number of gene filtering: ", dim(cm)[2])
  }
  ## Log2 transform
  message("Scaling and log2 transformation...")
  cm_norm = log(cm/scale_factor+1, base=log_base)
  if (filter_cell){
    ## Filter cells with low house keeping gene expressions 
    message("Filtering based on expressions of house keeping genes...")
    hg_list = hg_list[hg_list %in% rownames(cm)]
    hg_filtering = colMeans(cm_norm[hg_list,]) >= min_mean_hg
    cm_norm = cm_norm[,hg_filtering]
    cm = cm[,hg_filtering]
    message("Number of cells left after house-keeping gene filtering: ", dim(cm)[2])
  }
  ## Gene filtering
  message("Filtering gene with low expressions...")
  if (filter_gene_relaxed){
    genes_to_keep = rowSums(cm >= 2^gene_cutoff) >= gene_cell_cutoff
  }else{
    genes_to_keep = log2(rowMeans(cm)+1) >= gene_cutoff
  }
  cm = cm[genes_to_keep,]
  cm_norm = cm_norm[genes_to_keep,]
  message("The number of genes left after gene filtering: ", dim(cm)[1])
  ## Centering
  message("Centering gene expressions...")
  cm_norm_center = cm_norm-rowMeans(cm_norm)
  ## Check # of genes detected per cell
  if (verbose){
    message(paste0("Cell: ", dim(cm_norm)[2], "; ", "Gene: ", dim(cm_norm)[1]))
    message("Printing the number of genes expressed in all cells...")
    cat(summary(colSums(cm >= gene_detect)))
    cat("\n")
  }
  result = list()
  result[["raw_data"]] = cm
  result[["norm_data"]] = cm_norm
  result[["center_data"]] = cm_norm_center
  return (result)
}

## QC
## @param cm count matrix (raw or log transformed)
## @param hg_list list of house-keeping genes 
## @param scale_factor factor to scale tpm, default = 10
## @param log_base base of log, default = 2
## @param gene_detect detection limit for gene, default = 1
calculateQcMetrics <- function(cm, hg_list, scale_factor=10, log_base=2, gene_detect=1){
  ## Info for sample, plate and well
  sample = sapply(colnames(cm), function(x) unlist(strsplit(x, split=".", fixed=T))[1])
  plate = sapply(colnames(cm), function(x) unlist(strsplit(x, split=".", fixed=T))[2])
  well = sapply(colnames(cm), function(x) unlist(strsplit(x, split=".", fixed=T))[3])
  ## Total number of genes 
  gene = colSums(cm >= gene_detect)
  cm_norm = log2(cm/scale_factor+1)
  ## Average expression of house keeping genes
  hg_list = hg_list[hg_list %in% rownames(cm)]
  hk = colMeans(cm_norm[hg_list,])
  qc_df = cbind.data.frame(sample, plate, well, gene, hk)
  return(qc_df)
}

## Calculate pairwise cell correlation and tSNE embedding 
## @param cm_norm normalized count matrix
tSNEbyCor <- function(cm_norm, hc_method = "complete"){
  ## Compute pairwise pearson correlation and distance (1-r)
  message("Computing pairwise correlation and distance...")
  pairwise_cor = cor(cm_norm, method = "pearson")
  pairwise_dist = 1-pairwise_cor
  ## hirarchical clustering 
  message("Computing hierachical clustering...")
  hc = hclust(as.dist(pairwise_dist), method = hc_method)
  ## tSNE
  message("Computing tSNE embedding...")
  tsne_cor = Rtsne(as.dist(pairwise_dist), pca = F, is_distance=T)
  result = list()
  result[["pairwise_cor"]] = pairwise_cor
  result[["pairwise_dist"]] = pairwise_dist
  result[["hc_obj"]] = hc
  result[["tsne_obj"]] = tsne_cor
  return(result)
}

## Function to plot and store tSNE plots
## @param tsne_obj tSNE object 
## @param cm_norm log transformed cm
## @param color color scheme for tSNE
## @param title figure title
## @param filename file name of the figure 
## @param path path to store the figure
## @param dot_size size of dot on tSNE
## @param axis_size font size of axis labels
## @param legend_size font size of legend 
## @param save_obj whether to save the ggplot2 obj
## @param save_plot whether to save the plot
plotTsne <- function(tsne_obj, cm_norm, color, title="tSNE plot", filename="tSNE_cor.jpg", path="figures/",
                      dot_size=4, title_size=24, axis_size=18, legend_size=18, save_obj=FALSE, save_plot=TRUE){
  ## Check if length of tsne_obj and cm_norm matches
  try (if(dim(tsne_obj$Y)[1] != dim(cm_norm)[2]) stop("tsne_obj and cm_norm should have equal lengths..."))
   
  ## Generate a df for tSNE and color 
  tsne_mat = data.frame(tsne_obj$Y)
  rownames(tsne_mat) = colnames(cm_norm)
  colnames(tsne_mat) = c("tSNE1", "tSNE2")
  tsne_mat$color = color
  ## Generate ggplot2 obj
  g = ggplot(tsne_mat, aes_string(x= "tSNE1", y = "tSNE2", colour="color")) +
    geom_point(size = dot_size) + theme_bw() + ggtitle(title) + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = title_size),
          axis.title = element_text(size = axis_size), legend.text=element_text(size=legend_size))
  ## Display tSNE plot
  print(g)
  ## Save ggplot2 obj if save_obj = TRUE
  if (save_obj){
    return(g)
  }
  ## Save plot if save_plot = TRUE
  ggsave(filename, path = path, height = 8, width = 12)
}

## ==========================================================
## Functions for identifying highly variable genes
## ==========================================================

## Identify highly variable genes by pagoda2
## @param cm count matrix
## @param gam.k parameter for pagoda's adjustVariance, default = 10
## @param plot_var plot variance graph or not, default = TRUE
## @param n.cores # of cores to calculate variance, default = 6  
hvgPagoda <- function(cm, n.OdGenes=3000, gam.k=10, plot_var=TRUE, n.cores=6, verbose=TRUE){
  pagoda_obj <- Pagoda2$new(x = Matrix(as.matrix(cm), sparse=TRUE), n.cores=n.cores)
  ## Adjust the variance
  pagoda_obj$adjustVariance(plot=plot_var, gam.k=gam.k)
  ##pagoda_hvg = pagoda_obj$getOdGenes(n.odgenes=n.OdGenes)
  pagoda_hvg = pagoda_obj$misc$varinfo
  pagoda_hvg = pagoda_hvg[order(pagoda_hvg$lpa), ]
  if (verbose){
    cat(paste0("The number of highly variable genes identified by Pagoda2: ", length(pagoda_obj$misc$odgenes)))
  }
  return(pagoda_hvg)
}

## Identify highly variable genes by seurat 
## @param cm count matrix
## @param scale.factor default = 100000
## @param x.low.cutoff lower bound for expressions of highly variable genes, default = 0.0125
## @param x.high.cutoff upper bound for expressions of highly variable genes, default = 8
## @param y.cutoff lower bound for variance of highly variable genes, default = 1 
hvgSeurat <- function(cm, scale.factor=1E5, x.low.cutoff=0.0125, x.high.cutoff=8, y.cutoff=1, verbose=TRUE){
  seurat_obj = CreateSeuratObject(cm)
  ## Log normalize
  seurat_obj <- NormalizeData(
    object = seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = scale.factor
  )
  ## Detection of variable genes across the single cells
  seurat_obj <- FindVariableGenes(
    object = seurat_obj,
    mean.function = ExpMean, 
    dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, ## Default is 0.0125 
    x.high.cutoff = 8, ## Default is 8
    y.cutoff = 1 ## Default is 1
  )
  seurat_hvg = seurat_obj@var.genes
  cat(paste0("The number of highly variable genes identified by seurat: ", length(seurat_hvg)))
  return(seurat_hvg)
}

## ==========================================================
## Single cell signature scores  
## ==========================================================

## Compute plain average expression or corrected signature score 
## @param X.center, centered relative expression
## @param X.mean, average of relative expression of each gene (log2 transformed)
## @param n, number of genes with closest average expression for control genesets, default = 100 
## @param simple, whether use average without background correction, default  = FALSE
## @return return a vector of single cell score
scoreSignature <- function(X.center, X.mean, s, n=100, simple = FALSE, verbose=FALSE) {
  if(verbose) {
    message("cells: ", ncol(X.center))
    message("genes: ", nrow(X.center))
    message("genes in signature: ", length(s))
    message("Using simple average?", simple)
    message("processing...")
  }
  
  s <- intersect(rownames(X.center), s)
  message("genes in signature, and also in this dataset: ", length(s))
  ##message("These genes are: ", s)
  
  if (simple){
    s.score <- colMeans(X.center[s,])
  }else{
    s.score <- colMeans(do.call(rbind, lapply(s, function(g) {
      # g <- s[2]
      # message(g)
      if(verbose) message(".", appendLF = FALSE)
      g.n <- names(sort(abs(X.mean[g] - X.mean))[2:(n+1)])
      X.center[g, ] - colMeans(X.center[g.n, ])
    })))
  }
  
  if(verbose) message(" done")
  return(s.score)
}

## ==========================================================
## Functions for constructing seurat objects
## ==========================================================

## Construct a seurat object
## @para: count_matrix, cm
## @para: mito_symbol, mitocondrial genes in gene name/symbol format
## @para: min.cells, minimal number of cells that express a certain gene for that gene to be considered as expressed (default = 1)
## @para: is.expr, gene expression detection limit; default = 0 (for UMI); reads should use 5
## @para: project, project name for this seurat object 
seuratObj <- function(count_matrix, mito_symbol, min.cells = 1, is.expr = 0, project){
  seurat_obj <- CreateSeuratObject(round(count_matrix), min.cells = min.cells, is.expr = is.expr, project = project)
  mito_symbol = mito_symbol[mito_symbol %in% rownames(seurat_obj@raw.data)]
  percent.mito <- Matrix::colSums(seurat_obj@raw.data[mito_symbol, ]) / Matrix::colSums(seurat_obj@raw.data)
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.mito, col.name = "percent.mito")
}

## Create seurat obj, normalize, find variable genes, and scale data 
## @param cm, count matrix
## @param project project name
## @param min.cells minimal number of cells required to express a gene
## @param min.genes minimal number of genes required to be expressed in a cell
## @param is.expr detection limit for a gene
## @param scale.factor scale factor normalization, default=1E5
## @param do.scale whether to scale the data (divided by sd), default=F
## @param do.center whether to center the data (substract mean), default=T
preprocessSeuratObject <-function(cm, project, min.cells=0, min.genes=0, is.expr=0, 
                                  scale.factor=1E5, do.scale=F, do.center=T){
  ## Create Seurat Obj
  seurat_obj = CreateSeuratObject(cm, min.cells=min.cells, min.genes = min.genes, is.expr=is.expr, project=project)
  
  ## Normalize data
  seurat_obj <- NormalizeData(
    object = seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = scale.factor
  )
  
  ## Detection of variable genes across the single cells
  seurat_obj <- FindVariableGenes(
    object = seurat_obj,
    mean.function = ExpMean, 
    dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, ## Default is 0.0125 
    x.high.cutoff = 8, ## Default is 8
    y.cutoff = 1 ## Default is 1
  )
  ##length(x = seurat_obj@var.genes)
  
  ## Scaling the data and removing unwanted sources of variation
  seurat_obj <- ScaleData(
    object = seurat_obj,
    do.scale = do.scale,
    do.center = do.center
    #vars.to.regress = c("nUMI")
  )
  
  return(seurat_obj)
}

## Create seurat obj, normalize, find variable genes, and scale data 
## @param cm, count matrix
## @param project project name
## @param min.cells minimal number of cells required to express a gene
## @param min.genes minimal number of genes required to be expressed in a cell
## @param init_names_field, which field in the cell barcode name represents cell name 
## @param init_names_delim, deliminator to separate cell barcode name 
## @param scale.factor scale factor normalization, default=1E5
## @param do.scale whether to scale the data (divided by sd), default=F
## @param do.center whether to center the data (substract mean), default=T
preprocessSeuratObjectV3 <-function(cm, project, min.cells=0, min.genes=0, 
                                    init_names_field = 1, init_names_delim = ".",
                                    scale.factor=1E5, do.scale=F, do.center=T){
    ## Create Seurat Obj
    seurat_obj = CreateSeuratObject(cm, 
                                    min.cells=min.cells, min.features = min.genes, 
                                    project=project, 
                                    names.field = init_names_field, 
                                    names.delim = init_names_field)
    
    ## Normalize data
    seurat_obj <- NormalizeData(
        object = seurat_obj,
        normalization.method = "LogNormalize",
        scale.factor = scale.factor
    )
    
    ## Detection of variable genes across the single cells
    seurat_obj <- FindVariableFeatures(
        object = seurat_obj,
        selection.method = "vst"
    )
    ##length(x = seurat_obj@var.genes)
    
    ## Scaling the data and removing unwanted sources of variation
    all.genes <- rownames(cm)
    seurat_obj <- ScaleData(
        object = seurat_obj,
        features = all.genes,
        do.scale = do.scale,
        do.center = do.center
        #vars.to.regress = c("nUMI")
    )
    
    return(seurat_obj)
}

generatePseudobulk <- function(nmf_score, seurat_obj, out_path, aggr_method = "median", cutoff = 1){
    metagene_program_names = sort(unique(nmf_score$signature_1))
    pseudobulk = NULL
    for (metagene in metagene_program_names){
        cm = seurat_obj@raw.data[,nmf_score$signature_1 == metagene & nmf_score$score_1 >= 1]
        cm = as.matrix(cm)
        if (aggr_method == "mean"){
            pseudobulk = cbind(pseudobulk, rowMeans(cm))
        }
        if (aggr_method == "median"){
            pseudobulk = cbind(pseudobulk, rowMedians(cm))
        }
    }
    colnames(pseudobulk) = metagene_program_names
    rownames(pseudobulk) = rownames(seurat_obj@raw.data)
    saveRDS(pseudobulk, file=paste0(out_path, "pseudobulk_cm_", aggr_method, ".rds"))
}

addMitoProp <- function(seurat_obj, mito_symbol, V3 = T){
  if (V3){
    tmp = GetAssayData(seurat_obj, slot = "counts")
    mito_symbol = mito_symbol[mito_symbol %in% rownames(tmp)]
    percent.mito <- Matrix::colSums(tmp[mito_symbol, ]) /
      Matrix::colSums(tmp)
    seurat_obj <- AddMetaData(object = seurat_obj, 
                              metadata = percent.mito, 
                              col.name = "percent.mito")
  } else{
    mito_symbol = mito_symbol[mito_symbol %in% rownames(seurat_obj@raw.data)]
    percent.mito <- Matrix::colSums(seurat_obj@raw.data[mito_symbol, ]) / Matrix::colSums(seurat_obj@raw.data)
    seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.mito, col.name = "percent.mito")
  }
  return(seurat_obj)
}

## ==========================================================
## Misc functions 
## ==========================================================

## Preprocess ss2 data for pagoda analysis 
## @param cm count matrix
## @param total_gene_cutoff minimum total number of genes in a cell
## @param gene_expr gene detection limit for each cell 
## @param gene_per_cell minimum number of cells with expressions of a gene 
## @return filtered count matrix 
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

gene_symbol_to_ensembl_id <- function(gene, dataset){
  ensembl = useMart("ensembl")
  mart <- useDataset(dataset, useMart("ensembl"))
  ids <- getBM(filters="external_gene_name", attributes=c("ensembl_gene_id", "external_gene_name"), values=gene, mart=mart)
  non_duplicates <- which(duplicated(ids$external_gene_name) == FALSE)
  ids <- ids[non_duplicates, ] 
}

go_analysis <- function(sigOE_genes, allOE_genes, keyType="ENSEMBL", 
                        OrgDb="org.Hs.eg.db", ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05){
  ego <- enrichGO(gene = sigOE_genes, 
                  universe = allOE_genes, 
                  keyType = keyType, 
                  OrgDb = OrgDb, 
                  ont = ont, 
                  pAdjustMethod = pAdjustMethod, 
                  qvalueCutoff = qvalueCutoff, 
                  readable = TRUE)
  
  cluster_summary <- subset(data.frame(ego), select = -c(geneID))
  cluster_summary <- data.frame(ego)
  return(list(ego=ego, cluster_summary=cluster_summary))
}