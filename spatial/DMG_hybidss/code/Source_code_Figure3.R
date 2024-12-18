#################
### Figure 3a ###
#################
DimPlot(object = seurat_obj, 
        group.by = "annot_new", 
        cols = color_scheme,
        label = T, 
        pt.size = 3, 
        label.size = 0) & NoAxes()

#################
### Figure 3b ###
#################
genes_to_plot = c("PDGFRA", "SOX10", "OLIG1", "OLIG2")
Idents(seurat_obj) = seurat_obj$signature_1
tmp = VlnPlot(seurat_obj, 
              features = genes_to_plot,
              idents = c(paste0("OPC-like-", 1:3), "AC-like"),
              cols = color_scheme,
              group.by = "signature_1", 
              pt.size = 0.25, 
              sort = T, ncol = 4)
tmp & xlab("") & ylab("Log expression level") &
  theme(plot.title = element_text(size=16, face="plain"),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=12, angle=45, hjust=1),
        axis.text.y = element_text(size=12))

#################
### Figure 3c ###
#################
genes = c("ASCL1", "EGFR", "HES6", "BTG2", "DLL1", 
          "MFNG", "NNAT", "PCP4",
          "PDGFRA", "CSPG4", "GPR17", "EPN2", "PLLP")
genes = rev(genes)

## Compute pseudobulk 
Idents(seurat_obj) = seurat_obj$signature_1
pb = AverageExpression(seurat_obj, slot = "counts")
pb = pb$RNA
pb = log2(pb/10+1)
pb = pb - rowMeans(pb)

## Convert to df for heatmap 
df = pb[genes, paste0("OPC-like-", 3:1)]
df$gene = factor(rownames(df), level = genes)
df = melt(df)
colnames(df) = c("gene", "metaprogram", "expression")

## Plot heatmap
colors = rev(brewer.pal(9, "RdBu"))
ggplot(df, aes(x=metaprogram, y=gene, fill=expression)) + 
  geom_tile(color="black", size=1) + 
  scale_fill_gradientn(colors=colors) + 
  labs(x="", y="", fill="Relative\nExpression") +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=16, angle=45, hjust=1),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))

###################
### Figure 3d-f ###
###################
normal_OPCs = c("Pri-OPC", "OPC", "OAPC")
malignant_OPCs = metagene_order[1:3]

df = pairwise_cor_long
df = df[df$Cell_type %in% normal_OPCs & df$Metaprogram %in% malignant_OPCs,]

ggplot(df, aes(x=Metaprogram, y=Cell_type)) + 
  geom_point(aes(color=tumor_score, size=normal_score)) + scale_size_area(max_size=12) +
  scale_color_gradientn(colors=brewer.pal(n=9, name="Reds")) +
  labs(title = ref_name,
       x = "Malignant metaprogram",
       y = "Normal cell type",
       color = "Expression score\n(tumor cells)", 
       size = "Expression score\n(normal cells)") +
  theme_classic() + 
  theme(plot.title = element_text(size=24, hjust = 0.5),
        axis.title = element_text(size=24),
        axis.text.x=element_text(size=20, hjust=1, angle=45),
        axis.text.y=element_text(size=20),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed', colour = "grey"))

#################
### Figure 3g ###
#################
RSS_res = RSS_res[, c(6, 7, 3)]
colnames(RSS_res) = paste0("OPC_like_", 1:3)
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
### Figure 3g ###
#################
tmp = df[df$metaprogram %in% paste0("OPC-like-", 1:3) & df$type == "fresh", ]
plotProportionDot(tmp[tmp$location != "spinal",], 
                  multi_vars = F, 
                  dotplot = T,
                  x_var = "location", 
                  y_var = "proportion", 
                  fill_var = "location", 
                  color_scheme = colors_location,
                  d=1,
                  stat_color = "black",
                  fname = paste0("thalamus_vs_pons_fresh_OPC.pdf"),
                  out_dir = seurat_fig_folder)