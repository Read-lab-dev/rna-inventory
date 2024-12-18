library(infercnv)

## Compute score for oligodendrocytes or immune cells 
computeOligoImmuneScore <- function(cm, markerGenes, scoreName, geneAsRow = TRUE){
  ## Convert input cm into a data frame
  ## Traverse if geneAsRow == TRUE
  if (geneAsRow){
    t.df.tumor.data <- as.data.frame(as.matrix(t(cm)))
  }
  else{
    t.df.tumor.data <- as.data.frame(as.matrix(cm))
  }
  
  ## Subset gene expressions of only marker genes
  ## Compute average
  markers = intersect(colnames(t.df.tumor.data), markerGenes)
  df = dplyr::mutate((t.df.tumor.data), cell = rownames((t.df.tumor.data))) %>%
    dplyr::select(cell,markers)
  df = data.frame(df, score = rowMeans(df[,2:(length(markers) + 1)]))  
  
  ## Rename rows and columns
  rownames(df) = df$cell
  df = dplyr::select(df, -cell)
  colnames(df) = c(markers, scoreName)
  
  return(df)
}

## Add sample annotation 
generateSampleAnnotation <- function(oligo_df, immune_df, orig_annotations, cell_names, threshold = 4){
  ## Combine oligo and immune scores and annotations 
  df = cbind.data.frame(oligo_df$oligo_score, immune_df$immune_score, orig_annotations)
  rownames(df) = cell_names
  colnames(df) = c("oligo_score", "immune_score", "annotation")
  df$malig_status = ifelse((df$oligo_score > threshold | df$immune_score > threshold), "non-malignant", "malignant")
  df$patient = sapply(rownames(df), function(x) unlist(strsplit(x,split=".",fixed=TRUE))[1])
  df$malig_anno = ifelse(df$malig_status == 'malignant', 
                         paste("malignant", df$patient, sep="_"),
                         paste0(df$annotation,"(non-malignant)"))
  df = dplyr::select(df, malig_anno)
  return(df)
}

generateSampleAnnotationTenX <- function(oligo_df, immune_df, orig_annotations, cell_names, patient_id, threshold = 4){
  ## Combine oligo and immune scores and annotations 
  df = cbind.data.frame(oligo_df$oligo_score, immune_df$immune_score, orig_annotations)
  rownames(df) = cell_names
  colnames(df) = c("oligo_score", "immune_score", "annotation")
  df$malig_status = ifelse((df$oligo_score > threshold | df$immune_score > threshold), "non-malignant", "malignant")
  df$malig_anno = ifelse(df$malig_status == 'malignant', 
                         paste("malignant", patient_id, sep="_"),
                         paste0(df$annotation,"(non-malignant)"))
  df = dplyr::select(df, malig_anno)
  return(df)
}

## Need testing 
addMalignancyStatus <- function(observations, references, seurat_obj, k = 2, labels = c('tumor', 'normal'), hclust_method = 'ward.D', changeNames = TRUE){
  ## Tumor cell re-annotation
  hc <- hclust(dist(t(observations)), method = hclust_method)
  cell.annot.tumor <- as.character(factor(cutree(hc, k = k), labels = labels))
  names(cell.annot.tumor) <- colnames(observations)
  
  ## Normal cell re-annotation
  cell.annot.normal <- rep("normal", dim(references)[2])
  names(cell.annot.normal) <- colnames(references)
  
  ## Concat both
  cell.annot = c(cell.annot.tumor, cell.annot.normal)
  if (changeNames == TRUE){
    names(cell.annot) = gsub(".", "-", names(cell.annot), fixed = TRUE)
  }
  
  ## Sort cell.names and store tumor status into seurat obj
  cell.annot = cell.annot[seurat_obj@cell.names]
  print(table(names(cell.annot) == seurat_obj@cell.names))
  seurat_obj@meta.data$malignancy = cell.annot
  
  return(seurat_obj)
}

#This function identifies cells as malignant based on (1) their sum of squares of CNV signals across the genome, and (2) based
#on the correlation with the CNVs defined by exome-seq (approach from Filbin et. al., Science). 
#@param Ecnv_smoothed is the output of infer.cnv function, 
#@param malignant is a vector of malignant cell names as defined by gene expression.
#@param genes_cn is the output of match.exome.cnv function. 
#@param sample_ident is pre-defined.
#The output of this function is a vector of cell names identified as being malignant based on the specified cutoffs. 
#The function also generates a plot of CNV signal correlation vs CNV signal strength for each sample.
#If genes_cn is not supplied, the method will instead use correlation with the average CNV signal across all cells.

id.malignant = function(Ecnv_smoothed, malignant_exprs, sample_ident, plot_path = "figures/", genes_cn = 0, signal_cutoff = 0.02, correlation_cutoff = 0.2)
{
  correlation = Ecnv_smoothed[,-c(1,2,3,4)]
  #Retrieve sample IDs.
  Sample = sample_ident
  
  #remove data for which there is no exome data.
  if (genes_cn != 0)
  { 
    correlation = correlation[,Sample %in% colnames(genes_cn)]
  }
  
  unique_samples = unique(Sample)
  
  #Here, check if each unique sample contains at least 3 malignant cells. If there are less than 3, consider
  #all cells of the sample to be non-malignant.
  tally = vector()
  for (i in 1:length(unique_samples))
  {
    check = colnames(correlation)[Sample == unique_samples[i]]
    tally[i] = sum(check %in% malignant_exprs)
  }
  
  
  #The following vector is boolean, stores the locations of all malignant cells by cnv correlation.
  malignant = vector()
  #The following dataframe is for storing cnv_correlations
  cnv_correlations = data.frame(Cell = character(), Correlation = numeric())
  
  #if exome cnv is not supplied, calculate correlation with average inferred CNV signal.
  if(genes_cn == 0)
  {
    for (i in 1:length(unique_samples))
    {  
      temp = correlation[,Sample == unique_samples[i]]
      
      #The ifelse here guards against samples that do not have >2 malignant cells.
      #For those samples, the cnv_correlation of all cells is set to 0.
      if (tally[i] > 2)
      { 
        temp_malignant_exprs = temp[,colnames(temp) %in% malignant_exprs]
        genes_cn = rowMeans(temp_malignant_exprs)
        temp_cnv_correlations = apply(temp, 2, function(x) cor(x, genes_cn, method = "spearman"))
      }else{
        temp_cnv_correlations = rep(0, ncol(temp))
      }
      
      #There are a few cells for which temp_cnv_correlation becomes "NA" becomes the CNV values for that cell are all 0s.
      temp_cnv_correlations[is.na(temp_cnv_correlations)] = 0
      
      temp_cnv_correlations = data.frame(Cell = colnames(temp), Correlation = temp_cnv_correlations)
      cnv_correlations = rbind(cnv_correlations, temp_cnv_correlations)
      temp_malignant = temp_cnv_correlations$Correlation >= correlation_cutoff
      malignant = c(malignant,temp_malignant)
    }
    
    
  } else {
    for (i in 1:length(unique_samples))
    {
      #if exome cnv is supplied, calculate correlation with exome cnv as ground truth.
      temp = correlation[,Sample == unique_samples[i]]
      temp_exome = genes_cn[,colnames(genes_cn) == unique_samples[i]]
      temp_cnv_correlations = apply (temp, 2, function(x) cor(x, temp_exome, method = "spearman"))
      temp_cnv_correlations = data.frame(Cell = colnames(temp), Correlation = temp_cnv_correlations)
      cnv_correlations = rbind(cnv_correlations, temp_cnv_correlations)
      temp = temp_cnv_correlations$Correlation >= correlation_cutoff
      malignant = c(malignant, temp)
    }}
  
  #Calculate total cnv signal.
  cnv_signal = apply(correlation,2,function(x) sum(x^2)/length(x))
  
  #The following dataframe will be returned by the function.
  to_return = data.frame(Cell = colnames(correlation), CNV_signal = cnv_signal)
  to_return = inner_join(to_return, cnv_correlations, by = "Cell")
  
  #Filter by cnv_signal, and correlation.
  cnv_signal_malignant = to_return$CNV_signal >= signal_cutoff
  cnv_correlation_malignant = to_return$Correlation >= correlation_cutoff
  malignant = cnv_signal_malignant & cnv_correlation_malignant
  indeterminate = xor(cnv_signal_malignant, cnv_correlation_malignant)
  to_return = data.frame(to_return, Malignant = malignant, Indeterminate = indeterminate)
  
  #Generate plot of CNV correlation vs CNV signal for all cells.
  pdf(paste0(plot_path, "All samples CNV QC plot.pdf"), width = 7,height=7)
  toplot = data.frame(CNV_signal = to_return$CNV_signal,Correlation = to_return$Correlation)
  plot = ggplot(toplot,aes(x = CNV_signal, y = Correlation)) + geom_point()
  plot = plot + xlab("CNV signal") + ylab("CNV correlation")
  plot = plot + xlim(c(0, 0.2))
  plot = plot + geom_hline(yintercept = correlation_cutoff, color = "red")
  plot = plot + geom_vline(xintercept = signal_cutoff, color = "red")
  print(plot)
  dev.off()
  
  #Generate plot of CNV correlation vs CNV signal for all cells
  for (i in 1:length(unique_samples))
  {
    pdf(paste0(plot_path, unique_samples[i]," CNV QC plot.pdf"), width = 7,height=7)
    toplot = to_return[Sample == unique_samples[i],]
    toplot = data.frame(CNV_signal = toplot$CNV_signal, Correlation = toplot$Correlation)
    plot = ggplot(toplot,aes(x = CNV_signal, y = Correlation)) + geom_point()
    plot = plot + xlab("CNV signal") + ylab("CNV correlation")
    plot = plot + xlim(c(0, 0.2))
    plot = plot + geom_hline(yintercept = correlation_cutoff, color = "red")
    plot = plot + geom_vline(xintercept = signal_cutoff, color = "red")
    print(plot)
    dev.off()
    
  }
  
  return (to_return)
}
