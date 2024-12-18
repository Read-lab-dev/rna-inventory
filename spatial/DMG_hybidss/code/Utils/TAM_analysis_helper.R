#######################
### Geneset scoring ###
#######################

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

## Compute single cell gene signature score for a list of genesets  
## @param cm_center, centered relative expression
## @param cm_mean, average of relative expression of each gene (log2 transformed)
## @param nmf_gene_list, a list of genesets  
## @param simple, whether use average without background correction, default  = FALSE
## @return a matrix of single cell scores; row = geneset, column = cell id 
scoreNmfGenes <- function(cm_center, cm_mean, nmf_gene_list, verbose=TRUE, simple=FALSE){
  results = NULL
  row_names = NULL
  
  i = 1  
  for (genes_list in nmf_gene_list){
    nmf_name = names(nmf_gene_list)[i]
    row_names = c(row_names, nmf_name)
    scores = scoreSignature(cm_center, cm_mean, genes_list, verbose=verbose, simple=simple)
    results = rbind(results, scores)
    i = i + 1
  }
  
  colnames(results) = colnames(cm_center)
  rownames(results) = row_names
  ##results = results[sort(rownames(results)),]
  return(results)
}

################
### Plotting ###
################

## Plot proportions of cells grouped by a covariate as pie charts
## Based on a seurat object that store all the metadata 
## @param var_index, index of the condition to plot
## @param var, variable name of the covariate stored in the seurat object, default = sample (sample name)
## @param var_order, order of the conditions to plot, default = sample_order 
## @param nrow, the number of rows in the final plot
## @return print out or plot the resulting pie chart(s)
plotPieChart <- function(var_index, var="sample", var_order=sample_order, nrow=2){
  pie_chart_list = list()
  for (sample in var_order[var_index]){
    print(sample)
    to_subset = seurat_obj@meta.data[, var] == sample
    tmp2 = clust_order[clust_order %in% seurat_obj$seurat_clusters[to_subset]]
    crosstab = table(seurat_obj$seurat_clusters[to_subset])
    crosstab = crosstab / sum(crosstab)
    crosstab = data.frame(crosstab)
    colnames(crosstab) = c("Cluster", "Proportion")
    tmp = ggplot(crosstab, aes(x="", y=Proportion, fill=Cluster)) + 
      geom_bar(width = 2, color = "black", stat = "identity") + 
      coord_polar("y", start=0) + 
      labs(title=sample, x="", y="") + 
      scale_fill_manual(values = clust_colors) +
      theme(axis.line = element_blank(),
            axis.text = element_blank(),
            legend.position = "None")
    tmp
    pie_chart_list[[sample]] = tmp
  }
  plot_grid(plotlist = pie_chart_list, nrow = nrow)
}

## Plot -log10(pvalues) of over-represented GO terms by a hypergeometric test 
## Based on a list of GO term enrichment test reults by ClusterProfiler  
## @param target_cluster, name of the cell cluster to plot 
## @param go_term, which GO terms to show. This can be a regular expression 
## @param Add, whether to add the "myeloid cell activation involved in immune response" GO term, default=FALSE
## @param N, top N GO terms to show, default = 5
## @param nChar, the threshold of character numbers to wraps GO terms, default = 40
## @param width, @height, width and height of the final plot, default = 6 and 4
## @return NULL, save the final plots into a designated folder 
plotGOres <- function(target_cluster, go_term, Add=F, N=5, nChar=40, width=6, height=4){
  res = go_result[[paste0("c", target_cluster)]]$cluster_summary
  tmp = res[res$Description == "myeloid cell activation involved in immune response", c(1:7, 9)]
  df = res[grep(go_term, res$Description), c(1:7, 9)]
  if (Add){
    df = rbind(tmp, df)
  }
  df$log10p = -log10(df$p.adjust)
  
  df$Description = sapply(df$Description, function(x) {
    if(nchar(x)>nChar){
      tmp2 = unlist(stri_split_fixed(x, " "))
      len = length(tmp2)
      tmp3 = c(tmp2[1:floor(len/2)], "\n", tmp2[(floor(len/2)+1):len])
      tmp3 = paste(tmp3, collapse=" ")
      return(tmp3)
    } else{
      return(x)
    }
  })
  
  ggplot(df[1:N, ], aes(x=Description, y=log10p)) + 
    geom_col(color="black", fill=clust_colors[target_cluster]) + coord_flip() + 
    labs(x="", y="-log10(p value)") + 
    theme(axis.title = element_text(size=16),
          axis.text.x=element_text(size=16),
          axis.text.y=element_text(size=12),
          legend.position = "None")
  ggsave(paste0("pvalue_go_res_", target_cluster, ".pdf"), path=seurat_fig_folder, width=width, height=height)
}

## UMAP with cells of designated conditions shown 
## Based on a seurat object, in this case onlu focused on the treatment_merged metadata. Can be expanded to work with other metadata
## @param x, condition to plot 
## @param dot_size, size of the dots in the UMAP, default = 2
## @param alpha, transparency of the dots in the UMAP, default = 1
## @return NULL, save the final plots into a designated folder 
plotUMAPsub <- function(x, dot_size=2, alpha=1){
  Idents(seurat_obj) = seurat_obj$patient_id
  p = DimPlot(seurat_obj, 
              cells=colnames(seurat_obj)[seurat_obj$treatment_merged == x], 
              cols=patient_colors, pt.size=dot_size, label=F) & NoAxes() &  
    theme(legend.text = element_text(size=20))
  p[[1]]$layers[[1]]$aes_params$alpha = alpha
  p
  ggsave(paste0("Figure_6B_umap_", x, ".pdf"), path=seurat_fig_folder, 
         width=8, height=6)
}