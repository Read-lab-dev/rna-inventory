library(ggplot2)

## Make a bar graph to summarize proportion of clusters/programs/etc in each sample 
## @para: x, x-axis data
## @para: y, y-axis data
## @para: x_order: order of x-axis variables
## @para: y_order: order of y-axis variables
## @para: col_names: column names for df for plotting 
## @para: x_var: name of x-axis variable to plot 
## @para: y_var: name of y-axis variable to plot 
## @para: fill_var: name of variable for filling the bar graph 
## @para: plot_title: title of the plot
## @para: x_title: title of x-axis
## @para: y_title: title of y-axis
## @para: bar_position: stack vs dodge 
## The rest parameters are sizes of different labels  
plotProportion <- function(x, y, x_order, y_order, col_names, 
                           x_var, y_var, fill_var, colors,
                           plot_title="", x_title="", y_title="",
                           bar_position = "stack",
                           bar_size = 2,
                           axis_color = "black", axis_size = 2,
                           title_size=32, axis_title_size=28, 
                           x_text_size=20, x_text_angle=45, x_text_hjust=1,
                           y_text_size=24, legend_title_size=28, legend_text_size=24){
    crosstab = table(x, y)
    crosstab = crosstab/rowSums(crosstab)
    crosstab = crosstab[x_order, y_order]
    crosstab = data.frame(crosstab)
    print(crosstab)
    
    colnames(crosstab) = col_names
    ggplot(crosstab, aes_string(x=x_var, y=y_var, fill=fill_var)) + 
        geom_bar(stat="identity", color="black", size = bar_size,
                 position = bar_position) + 
        scale_fill_manual(values = colors) + 
        ggtitle(plot_title) + xlab(x_title) + ylab(y_title) + 
        theme_classic() + 
      theme(axis.line.x = element_line(color = axis_color, size = axis_size),
            axis.line.y = element_line(color = axis_color, size = axis_size),
            plot.title = element_text(hjust=0.5, face="bold", size=title_size), 
            axis.title = element_text(size=axis_title_size),
            axis.text.x = element_text(angle=x_text_angle, hjust=x_text_hjust, size=x_text_size),
            axis.text.y = element_text(size=y_text_size),
            legend.title = element_text(size=legend_title_size),
            legend.text = element_text(size=legend_text_size))
}

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

## Function to plot cell types with subclusters
plotSubcluster <- function(cell_type, annot_col, color_scheme = clust_colors){
  col_name = paste0('c', cell_type, "_sub")
  seurat_obj@meta.data[,col_name] = sapply(seurat_obj@meta.data[,annot_col], function(x){
    x = as.character(x)
    ifelse(startsWith(x, cell_type), x, "Others")
  })
  
  colors = c(color_scheme[grep(cell_type, names(color_scheme))],
             "Others"  = "lightgrey")
  ##names(colors)[1] = col_name
  
  DimPlot(seurat_obj, 
          group.by = col_name, 
          cols = colors,
          label = TRUE, 
          pt.size = 1, 
          label.size = 0) + 
    NoAxes() + ggtitle("") + 
    theme(legend.text = element_text(size=20), legend.position = c(0.1, 0.1))
  ggsave(paste0("UMAP_", col_name, ".pdf"), 
         path=seurat_fig_folder, width=8, height=8)
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

## Make a grid of violin plots for gene expressions for each cell subpopulation 
## @para: cm, count matrix (log transformed value)
## @para: annotation, cell annotation of each cell 
## @para: genes, genes to plot
## @para: drug, name of the drug (only used for image name)
## @para: out, output directory
## The rest label sizes and image sizes 
plotGridViolin <- function(cm, annotation, genes, drug, out = seurat_fig_folder,
                           y_axis_title = 24, x_axis_text = 24, strip_text = 12,
                           out_width = 16, out_height = 12){
  tmp = data.frame(cm[,genes])
  tmp$metaprogram = annotation
  tmp = melt(tmp, id.vars = "metaprogram")
  colnames(tmp) = c("metaprogram", "gene", "expression")
  
  ggplot(tmp, aes(x=metaprogram, y=expression, fill = metaprogram)) + 
    geom_violin(scale = "width") + 
    xlab("") + ylab("Log2 expression\n") + 
    scale_fill_manual(values = color_scheme) +
    facet_grid(gene~.) + 
    theme(axis.title.y = element_text(size = y_axis_title, face = "bold"),
          axis.text.x = element_text(size = x_axis_text, face = "bold", angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          strip.text = element_text(size=strip_text, face="bold"),
          axis.ticks.y = element_blank(),
          plot.margin = margin(10,10,10,70),
          legend.position = "none")
  ggsave(paste0("PF_", drug, "_targets_violin.png"), 
         path = out, width = out_width, height = out_height)
}

## Plot dotplot with p value and jaccard index
plotProgramOverlap <- function(a, b, a_name, b_name, x_lab, y_lab,
                               color_lab = "-log10(p-value)", size_lab = "Jaccard Index",
                               title_size=32, axis_title_size=28, 
                               x_text_size=20, x_text_angle=45, x_text_hjust=1,
                               y_text_size=24, legend_title_size=28, legend_text_size=24){
    ## Compute pairwise fisher exact pvalue and jaccard index
    bg_genes = unique(c(unlist(a), unlist(b)))
    pvalue_res = NULL
    jaccard_res = NULL
    
    for (prog1 in a){
        tmp = NULL
        tmp2 = NULL
        for (prog2 in b){
            tmp = c(tmp, fisher_test(prog1, prog2, bg_genes))
            tmp2 = c(tmp2, jaccard_index(prog1, prog2))
        }
        pvalue_res = rbind(pvalue_res, tmp)
        jaccard_res = rbind(jaccard_res, tmp2)
    }
    
    rownames(pvalue_res) = names(a); colnames(pvalue_res) = names(b)
    rownames(jaccard_res) = names(a); colnames(jaccard_res) = names(b)
    
    ## Combine p value and jaccard index in long format 
    pvalue_res = melt(pvalue_res)
    jaccard_res = melt(jaccard_res)
    df = cbind.data.frame(pvalue_res$Var1, pvalue_res$Var2, pvalue_res$value, jaccard_res$value)
    colnames(df) = c(a_name, b_name, "P_value", "Jaccard_Index")
    df$P_value = -log10(df$P_value)
    
    ## Plot dotplot
    ggplot(df, aes_string(x=a_name, y=b_name)) +
        geom_point(aes(color=P_value, size=Jaccard_Index)) + 
        scale_size_area(max_size=20) +
        scale_color_gradientn(colors=(brewer.pal(n=9, name="YlOrRd"))) +
        xlab(x_lab) + ylab(y_lab) +
        labs(color=color_lab, size=size_lab) +
        theme_classic() + 
        theme(plot.title = element_text(hjust=0.5, face="bold", size=title_size), 
              axis.title = element_text(face="bold", size=axis_title_size),
              axis.text.x = element_text(angle=x_text_angle, hjust=x_text_hjust, size=x_text_size),
              axis.text.y = element_text(size=y_text_size),
              legend.title = element_text(face="bold", size=legend_title_size),
              legend.text = element_text(size=legend_text_size),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed', colour = "grey"),
              plot.margin = margin(10,10,10,10))
}

## Heatmap for metaprogram scores 
## @para: nmf_score, metaprogram score (row = cells, col = scores of each metaprogram)
## @para: annotation, metaprogram annotation of each cell 
## @para: metagene_order, order of metaprograms to plot 
## @para: cutoff, upper and lower cutoff to trim values, default = 2.5
## @para: col_pal, color pallete to use, default = "RdBu" 
## @para: col_range, range of color to use from the color pallete, default = 1:9
## @para: out_dir, name of the output directory
## @para: out_name, name of the output figure file, default = "nmf_cell_score_sorted.png"
## @para: font_size, font size of labels of heatmap, default = 16
plotMetaScore <- function(nmf_score, annotation, metagene_order, 
                          cutoff=2.5, col_pal="RdBu", col_range=1:9,
                          annot_color_row = color_scheme, 
                          annot_color_col = color_scheme, 
                          out_dir = nmf_fig_folder,
                          out_name = "nmf_cell_score_sorted",
                          plot_pdf = F,
                          width=12, height=8){
    
  ## Order cells based on assigned metaprograms
  cell_order = annotation
  names(cell_order) = rownames(nmf_score)
  cell_order = factor(cell_order, levels = metagene_order)
  cell_order = cell_order[order(cell_order)] 
  cell_order_df = data.frame(cell_order)
  colnames(cell_order_df) = "annotation"
  
  ## Side bar for rows
  metagene_order_df = data.frame(metagene_order)
  rownames(metagene_order_df) = metagene_order
  colnames(metagene_order_df) = "metaprogram"
    
  ## Sort cells and metaprograms (cells: expressed program; metaprograms: specified order)
  nmf_score_sorted = nmf_score[names(cell_order), metagene_order]
  nmf_score_sorted = as.matrix(nmf_score_sorted)
  
  ## Trim values
  nmf_score_sorted_heatmap = ifelse(nmf_score_sorted > cutoff, cutoff, 
                                    ifelse(nmf_score_sorted < -cutoff, -cutoff, 
                                           nmf_score_sorted))
  
  ## heatmap for NMF scores
  hm_colors = rev(brewer.pal(n=9, name=col_pal))[col_range]
  hm_colors = colorRampPalette(colors = hm_colors)(100)
  fname = paste0(out_dir, out_name)
  if (plot_pdf){
    pdf(file=paste0(fname, ".pdf"), width = width, height = height)
  } else{
    png(filename=paste0(fname, ".png"), width = width*100, height = height*100)
  }
  pheatmap(t(nmf_score_sorted_heatmap), 
           color = hm_colors, 
           cluster_rows = F, cluster_cols = F, 
           labels_col = "", labels_row = "",
           annotation_row = metagene_order_df,
           annotation_col = cell_order_df, 
           annotation_colors = list(annotation = annot_color_col, metaprogram=annot_color_row),
           annotation_names_row = F, 
           annotation_names_col = F)
  dev.off()
}

## Heatmap for metaprogram scores 
## @para: cm_center, centered count matrix (row = genes, col = cells)
## @para: nmf_score, metaprogram score (row = cells, col = scores of each metaprogram)
## @para: annotation, metaprogram annotation of each cell 
## @para: metagene_order, order of metaprograms to plot 
## @para: cutoff, upper and lower cutoff to trim values, default = 2.5
## @para: col_pal, color pallete to use, default = "RdBu" 
## @para: col_range, range of color to use from the color pallete, default = 1:9
## @para: out_dir, name of the output directory
## @para: out_name, name of the output figure file, default = "nmf_cell_score_sorted.png"
## @para: font_size, font size of labels of heatmap, default = 16
plotMetaGeneExpr <- function(cm_center, nmf_score, annotation, metagene_order, 
                             nmf_genes = nmf_marker_genes_final,
                             cutoff=2.5, col_pal="RdBu", col_range=1:9,
                             annot_color_row = color_scheme, 
                             annot_color_col = color_scheme, 
                             out_dir = nmf_fig_folder,
                             out_name = "nmf_marker_gene_expr",
                             plot_pdf = F,
                             width=12, height=8){
  ## Order cells based on assigned metaprograms
  cell_order = annotation
  names(cell_order) = rownames(nmf_score)
  cell_order = factor(cell_order, levels = metagene_order)
  cell_order = cell_order[order(cell_order)] 
  cell_order_df = data.frame(cell_order)
  colnames(cell_order_df) = "annotation"
  
  ## heatmap for expressions of genes in merged NMF factors 
  hm_colors = rev(brewer.pal(n=9, name=col_pal))[col_range]
  hm_colors = colorRampPalette(colors = hm_colors)(100)
  cm_center = as.matrix(cm_center)
  cm_heatmap = ifelse(cm_center > cutoff, cutoff, 
                      ifelse(cm_center < -cutoff, -cutoff, cm_center))
  
  ## Only keep genes that are also in this dataset 
  nmf_genes = lapply(nmf_genes, function(x) x[x %in% rownames(cm_center)])
  nmf_genes = nmf_genes[metagene_order]
  nmf_marker_genes = unlist(nmf_genes)
  
  ## Prepare row annotations 
  nmf_marker_gene_vec = lapply(names(nmf_genes), function(x) rep(x, length(nmf_genes[[x]])))
  nmf_marker_gene_vec = unlist(nmf_marker_gene_vec)
  nmf_marker_gene_df = data.frame(nmf_marker_gene_vec)
  colnames(nmf_marker_gene_df) = "gene"
  rownames(nmf_marker_gene_df) = 1:nrow(nmf_marker_gene_df)
  
  ## Subset and rename counts
  tmp = cm_heatmap[nmf_marker_genes, names(cell_order)]
  rownames(tmp) = 1:nrow(tmp)
  
  ## Heatmap
  fname = paste0(out_dir, out_name)
  if (plot_pdf){
    pdf(file=paste0(fname, ".pdf"), width = width, height = height)
  } else{
    png(filename=paste0(fname, ".png"), width = width*100, height = height*100)
  }
  pheatmap(tmp, 
           color = hm_colors, 
           cluster_rows = F, cluster_cols = F, 
           labels_col = "", labels_row = "",
           annotation_row = nmf_marker_gene_df,
           annotation_col = cell_order_df, 
           annotation_colors = list(annotation = annot_color_col, gene=annot_color_row),
           annotation_names_row = F, 
           annotation_names_col = F)
  dev.off()
}