library(NMF)

## preprocess cm for nmf analysis
## @param cm log transformed and centered cm
nmf_df_preprocessing <- function(cm){
  ## convert negative values to zero
  cm = ifelse(cm < 0, 0, cm)
  ## remove genes with zeros in all cells 
  cm = cm[Matrix::rowSums(cm) != 0,]
}

## Find genes with high NMF score 
## @param w w-matrix (gene x NMF-factor)
## @param sample the sample name for current sample 
## @param num_genes number of genes to return
## Consider convert this into a vector so that this can be used in all places where scores are computed for NMF genes 
findNmfCorrelatedGenes <- function(w, sample, return_list = T, num_genes=30, type="nmf"){
  ## Find top genes with highest nmf score
  ## Store in a list and return that list 
  message(paste0("Identify top ", num_genes, 
                 " genes with highest NMF scores for sample ", sample))
  nmf_gene_list = list()
  nmf_gene_vec = NULL
  
  for (i in seq(dim(w)[2])){
    genes = names(sort(w[,i], decreasing = T))[1:num_genes]
    tmp = paste0(sample, "_", type, "_", i)
    if (return_list){
      ##nmf_gene_list[[tmp]] = names(sort(nmf_gene_cor[,i], decreasing = T))[1:num_genes]
      nmf_gene_list[[tmp]] = genes 
    } else{
      tmp = paste(rep(tmp), seq(num_genes), sep = "_")
      names(genes) = tmp
      nmf_gene_vec = c(nmf_gene_vec, genes)
    }
  }
  
  if (return_list){
    return(nmf_gene_list)
  } else{
    return(nmf_gene_vec)
  }
}

concatNmfFactor <- function(nmf_obj_list, all_genes, nmf=TRUE){
  sample_name = names(nmf_obj_list)
  
  result = NULL
  col_names = NULL
  i = 1 
  for (nmf_obj in nmf_obj_list){
    if (nmf){
      current_basis = basis(nmf_obj)
    }else{
      current_basis = nmf_obj$theta
    }
    current_basis = as.matrix(current_basis)
    current_basis = addZeroGeneCounts(current_basis, all_genes)
    for (j in seq(dim(current_basis)[2])){
      nmf_factor_name = paste(sample_name[i], "nmf", j, sep = "_")
      col_names = c(col_names, nmf_factor_name)
    }
    ## cbind won't work for NULL and data.frame
    result = cbind(result, current_basis)
    i = i+1
  }
  colnames(result) = col_names
  return(result)
}

addZeroGeneCounts <- function(bm, all_genes){
  genes = rownames(bm)
  to_add_genes = setdiff(all_genes, genes)
  ##print(to_add_genes)
  nrow = length(to_add_genes)
  ncol = dim(bm)[2]
  to_add_m = matrix(rep(0, nrow*ncol), nrow=nrow, ncol=ncol)
  rownames(to_add_m) = to_add_genes 
  bm = rbind(bm, to_add_m)
  bm = bm[all_genes,]
  return(bm)
}

clusterNmfFactors <- function(nmf_basis, cor_method="pearson"){
  nmf_dist = 1-cor(nmf_basis, method=cor_method)
  hc_nmf = hclust(as.dist(nmf_dist), method="ward.D2")
  result = list()
  result[["cor_coef"]] = 1-nmf_dist
  result[["dist"]] = nmf_dist
  result[["hc_obj"]] = hc_nmf
  return(result)
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

computeAverageCor <- function(factors, mat){
    mat_subset = as.vector(mat[factors, factors])
    mat_subset = as.vector(mat_subset)
    return(mean(mat_subset[mat_subset != 1]))
}

mergeNmfGenes <- function(nmf_basis, nmf_factor_sets, n = 30){
  nmf_basis_subset = nmf_basis[ ,nmf_factor_sets]
  if (length(nmf_factor_sets) == 1){
    return(names(sort(nmf_basis_subset, decreasing=T))[1:n])
  }else{
    avg_gene = rowMeans(nmf_basis_subset)
    return(names(sort(avg_gene, decreasing=T))[1:n])
  }
}

calculateNmfProgramStats <- function(nmf_basis, nmf_factor_sets){
  nmf_basis_subset = nmf_basis[ ,nmf_factor_sets]
  if (length(nmf_factor_sets) == 1){
    return(sort(nmf_basis_subset, decreasing = T))
  }else{
    avg_gene = rowMeans(nmf_basis_subset)
    sd_gene = apply(nmf_basis_subset, 1, sd)
    cv = sd_gene/avg_gene
    res = cbind.data.frame(avg_gene, sd_gene, cv)
    colnames(res) = c("mean", "sd", "cv")
    return(res[order(res$mean, decreasing=T),])
  }
}

metagene_score_signature <- function(df, num_metagene, rank=1){
  col_names = colnames(df)
  df$tmp1 = apply(df[,1:num_metagene], 1, 
                  function(x) sort(x, decreasing = T)[rank])
  df$tmp2 = apply(df[,1:num_metagene], 1, 
                  function(x) names(sort(x, decreasing = T))[rank])
  col_names = c(col_names, paste0("score_", rank), paste0("signature_", rank))
  colnames(df) = col_names
  return(df)
}

multi_expr <- function(nmf_score){
  nmf_score$diff = nmf_score$score_1 - nmf_score$score_2
  nmf_score$coexpr = ifelse(nmf_score$score_1 > 1 & nmf_score$diff > 1, "unique",
                            ifelse(nmf_score$score_2 > 1 & nmf_score$diff < 1, "double", "others"))
  nmf_score$expr = ifelse(nmf_score$score_1 >= 1, T, F)
  return(nmf_score)
}

# findNmfCorrelatedGenes_obsolete <- function(cm, h, sample, use_cor = F, transposed=F, 
#                                    cor_method="spearman", num_genes=30){
#   ## Transpose count and h matrix if the input are not transposed 
#   if (!transposed){
#     message("Transpose cm and h matrix if cells are not rows...")
#     h_t = t(h)
#     cm_t = t(cm)
#   }
#   
#   ## Check if cm_t and h_t have equal number of rows 
#   try(if(dim(cm_t)[1] != dim(h_t)[1])
#     stop("count and h matrix should have equal number of rows..."))
#   
#   ## compute pairwise correlation between each gene and each nmf factor 
#   message("Compute pairwise correlation between genes and NMF factors...")
#   nmf_gene_cor = cor(cm_t, h_t, method = cor_method)
#   
#   ## Find top genes with highest nmf score
#   ## Store in a list and return that list 
#   message(paste0("Identify top ", num_genes, 
#                  " genes with highest NMF scores..."))
#   nmf_gene_list = list()
#   for (i in seq(dim(h_t)[2])){
#     tmp = paste0("nmf_", i)
#     nmf_gene_list[[tmp]] = names(sort(nmf_gene_cor[,i], decreasing = T))[1:num_genes]
#   }
#   
#   return(nmf_gene_list)
# }