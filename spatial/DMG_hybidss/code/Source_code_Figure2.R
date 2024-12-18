#################
### Figure 2a ###
#################
DimPlot(object = seurat_obj, 
        group.by = "annot_new", 
        cols = color_scheme,
        label = T, 
        pt.size = 3, 
        label.size = 0) & NoAxes()

#################
### Figure 2b ###
#################
genes_to_plot = c("TOP2A", "UBE2C", "MCM2", "MCM5",
                  "PDGFRA", "OLIG1", "EGR1", "ASCL1",
                  "BCAS1", "MBP", "TF", "MOG",
                  "APOE", "AQP4", "ALDOC", "GFAP",
                  "GAP43", "DDIT3", "TIMP1", "SPP1")
genes_to_plot = rev(genes_to_plot)

tmp = DotPlot(seurat_obj, 
              features = genes_to_plot,
              cols = "RdBu",
              group.by = "annot_new")
tmp + coord_flip() + 
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(size=20, angle=45, hjust=1),
        axis.text.y = element_text(size=20))

#################
### Figure 2c ###
#################
#' Heatmap for metaprogram scores 
#' @param cm_center, centered count matrix (row = genes, col = cells)
#' @param nmf_score, metaprogram score (row = cells, col = scores of each metaprogram)
#' @param annotation, metaprogram annotation of each cell 
#' @param metagene_order, order of metaprograms to plot 
#' @param cutoff, upper and lower cutoff to trim values, default = 2.5
#' @param col_pal, color pallete to use, default = "RdBu" 
#' @param col_range, range of color to use from the color pallete, default = 19
#' @param out_dir, name of the output directory
#' @param out_name, name of the output figure file, default = "nmf_cell_score_sorted.png"
#' @param font_size, font size of labels of heatmap, default = 16
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

nmf_score = nmf_score_final_t[, 1:8]
colnames(nmf_score) = names(nmf_marker_genes_final) 
annotations = nmf_score_final_t$signature_1 
annotations = plyr::mapvalues(annotations, old_names, new_names)
plotMetaGeneExpr(cm_center, nmf_score, annotations, metagene_order)

#################
### Figure 2d ###
#################
#' Make a bar graph to summarize proportion of clusters/programs/etc in each sample 
#' @param x, x-axis data
#' @param y, y-axis data
#' @param x_order order of x-axis variables
#' @param y_order order of y-axis variables
#' @param col_names column names for df for plotting 
#' @param x_var name of x-axis variable to plot 
#' @param y_var name of y-axis variable to plot 
#' @param fill_var name of variable for filling the bar graph 
#' @param plot_title title of the plot
#' @param x_title title of x-axis
#' @param y_title title of y-axis
#' @param bar_position stack vs dodge 
#' The rest parammeters are sizes of different labels  
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

sample_order = c("MUV31", "MUV32", "MUV91",
                 "MUV1", "MUV10", "MUV83",
                 "MUV5",  "MUV17", "MUV35", "MUV77",  
                 "MUV78", "MUV86", "MUV87", 
                 "BT1126", "BT836", "BT869", "BT1462", "MUV54")
annot_order = c("Cycling", "OPC-like", "OC-like", "AC-like", "MES-like")

tmp = plotProportion(seurat_obj@meta.data$sample, 
                     seurat_obj@meta.data$annot_new,
                     sample_order, annot_order,
                     c("Sample", "Metaprogram", "Proportion"),
                     "Sample", "Proportion", "Metaprogram", 
                     color_scheme,
                     y_title = "Proportion", 
                     axis_title_size = 24,
                     x_text_size = 24, y_text_size = 24,
                     legend_title_size = 24,
                     bar_size = 1.5, axis_size = 1.5)
tmp

#################
### Figure 2e ###
#################
# Plot a sorted point plot for each regulon within a cell type based on their RSS
#' @param metaprog metaprograms to plot
#' @param df RSS result df
#' @param title_size, size of the title
#' @param axis_title_size, size of axis titles
#' @param axis_text_size, size of axis text
#' @param label_size, size of the label for highlighted regulons 
plotRankedRegulon <- function(metaprog, df = RSS_res,
                              num_highlight = 5, 
                              title_size=36, axis_title_size=28,
                              axis_text_size=24, label_size=5){
  ## Subset target regulon and sort based on RSS
  tmp = data.frame(sort(df[,metaprog], decreasing = T))
  colnames(tmp) = metaprog
  ## Add order, whether to highlight (top num_highlight), name of highlighted regulons 
  tmp$order = 1:nrow(tmp)
  tmp$highlight = c(rep("Yes", num_highlight), 
                    rep("No", nrow(tmp) - num_highlight))
  tmp$regulon = sapply(1:nrow(tmp), 
                       function(x){
                         if (tmp$highlight[x] == "Yes"){
                           extractTFs(rownames(tmp)[x])
                         }else{
                           ""
                         }
                       })
  
  ## Plot
  ggplot(tmp, aes_string(x="order", y=metaprog, 
                         color="highlight", label="regulon")) + 
    geom_point(size=1.5) + geom_text_repel(size=label_size, segment.alpha=0.5) +
    ggtitle(gsub("_","-",metaprog)) + xlab("Regulon") + 
    ylab("Specificity score") + 
    scale_color_manual(values = c("Yes" = "red", "No" = "blue")) + 
    theme(plot.title = element_text(size=title_size),
          axis.title = element_text(size=axis_title_size),
          axis.text = element_text(size=axis_text_size),
          axis.text.x = element_text(size=0),
          legend.position = "None")
}

colnames(RSS_res) = gsub("-", "_", colnames(RSS_res), fixed = T)
num_cols = ncol(RSS_res)

## Plot ranked RSS scores for each cell type 
ranked_plots = lapply(colnames(RSS_res), 
                      function(x) plotRankedRegulon(x, num_highlight = 10,
                                                    title_size = 28,
                                                    axis_title_size = 20,
                                                    axis_text_size = 20,
                                                    label_size = 5))
plot_grid(plotlist = ranked_plots, ncol = num_cols)

#################
### Figure 2f ###
#################
plotProportionDot <- function(df, 
                              multi_vars = T,
                              dotplot = T,
                              cond_var = "location", 
                              cond_var_value = "thalamus",
                              x_var = "age", 
                              y_var = "proportion", 
                              fill_var = "age", 
                              color_scheme = colors_age,
                              d=2,
                              stat_color="red",
                              fname = paste0("ped_vs_adult_thalamus.png"),
                              out_dir = seurat_fig_folder){
  if (multi_vars){
    cond = df[, cond_var] == cond_var_value
  } else{
    cond = rep(T, ncol(df))
  }
  p = ggplot(df[cond,], aes_string(x=x_var, y=y_var, fill=x_var)) 
  if (dotplot){
    p = p + geom_dotplot(binaxis='y', stackdir='center', dotsize = d) + 
      stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                   geom="pointrange", color=stat_color) 
  }
  else{
    p = p + geom_boxplot()
  }
  p = p + scale_fill_manual(values = color_scheme) + 
    xlab("") +
    theme(strip.text = element_text(size=20),
          strip.background = element_blank(),
          axis.title = element_text(size=20),
          axis.text.x = element_text(size=20, hjust=1, angle=45),
          axis.text.y = element_text(size=20),
          legend.position = "None")
  p + facet_wrap(.~metaprogram, ncol = 4)
  ggsave(fname, path=out_dir, width=12, height=8)                     
}

plotProportionDot(df, 
                  multi_vars = F,
                  dotplot = plotDotPlot,
                  x_var = "age", 
                  y_var = "proportion", 
                  fill_var = "age", 
                  color_scheme = colors_age,
                  d=1.5,
                  stat_color = "black",
                  fname = paste0("ped_vs_adult_wo_DIPG67.pdf"),
                  out_dir = seurat_fig_folder)

#################
### Figure 2g ###
#################
plotProportionDot(df[df$location != "spinal",], 
                  multi_vars = F, 
                  dotplot = F,
                  x_var = "location", 
                  y_var = "proportion", 
                  fill_var = "location", 
                  color_scheme = colors_location,
                  d=1,
                  stat_color = "black",
                  fname = paste0("thalamus_vs_pons_wo_DIPG67.pdf"),
                  out_dir = seurat_fig_folder)

#################
### Figure 2i ###
#################
OPC_score = pmax(seurat_obj$OPC_like, seurat_obj$Ribo_active, seurat_obj$IER)
mature_score = pmax(seurat_obj$AC_like, seurat_obj$OC_like, seurat_obj$MES_like)
seurat_obj = AddMetaData(seurat_obj, OPC_score - mature_score, "stemness_score")
seurat_obj = AddMetaData(seurat_obj, seurat_obj$OC_like - seurat_obj$AC_like, "lineage_score")

ggplot(seurat_obj@meta.data, aes(x=lineage_score, y=stemness_score, color=metaprogram_new_no_cc)) + geom_point(alpha=0.5, size=2) +
  scale_color_manual(values=color_scheme) +
  labs(x="Lineage score", y="Stemness score", color="Metaprogram") +
  facet_wrap(.~age, ncol = 2) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=16),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        strip.background=element_blank(),
        strip.text=element_text(size=16))