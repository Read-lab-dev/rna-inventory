#################
### Figure 4a ###
#################
DimPlot(
  object = seurat_obj,
  group.by = 'annotation_sep',
  cols = color_scheme,
  pt.size = 2,
  label = F,
  repel = TRUE) + NoAxes() + 
  theme(legend.text = element_text(size=20))

#################
### Figure 4b ###
#################
genes = c("HSPB8", "KCNN3", "AQP4", "GFAP",
          "AKAP12", "GAP43", "CHI3L1", "ANXA2",
          "BCAS1", "SIRT2", "SOX10", "MBP", 
          "OLIG1", "OLIG2",
          "MYT1", "EPN2", "ASCL1", "CSPG4",
          "GABBR2", "GRIA1", "CAMK2B", "SPARCL1"
)
genes = rev(genes)

tmp = DotPlot(seurat_obj, 
              group.by = "annotation_sep", 
              features = genes, 
              cols = "RdBu")
tmp + coord_flip() + 
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(size=20, angle=45, hjust=1),
        axis.text.y = element_text(size=20))

#################
### Figure 4c ###
#################
tmp = cor_res_expr_peak %>% filter(r >= 0.2 & pvalue_adj <= 0.05)

## Parameters
annot_order = c("OPC-like", "OC-like", "AC-like", "MES-like", "Synaptic")
color_scheme = c(color_scheme, "Synaptic" = "orangered")
color_scheme = color_scheme[annot_order]

cutoff=3; col_pal="RdBu"; col_range=1:5
out_dir = seurat_fig_folder
out_name = "Gene_activity_heatmap"
out_name2 = "Gene_expr_heatmap"

## Prepare counts for heatmap
set.seed(123456)
pb_atac_z = pb_atac_norm[rownames(tmp), ]
pb_atac_z = t(scale(t(pb_atac_z)))
pb_rna_z = pb_rna_norm[tmp$gene, ]
pb_rna_z = t(scale(t(pb_rna_z)))
to_subset = sample(rownames(pb_atac_z), 5000)
pb_atac_z = pb_atac_z[to_subset, ]
pb_rna_z = pb_rna_z[tmp[to_subset, "gene"],]

## Trim values
pb_atac_z = ifelse(pb_atac_z > cutoff, cutoff, 
                   ifelse(pb_atac_z < -cutoff, -cutoff, pb_atac_z))
pb_rna_z = ifelse(pb_rna_z > cutoff, cutoff, 
                  ifelse(pb_rna_z < -cutoff, -cutoff, pb_rna_z))

## hclust genomic regions 
hc_res1 = hclust(dist(pb_atac_z), method = "ward.D2")
hc_res2 = hclust(dist(t(pb_atac_z)), method = "ward.D2")

## Cell order
##cell_order = meta[colnames(pb_atac_z), "annotation_sep"]
cell_order = cells_NNs_annot
names(cell_order) = colnames(pb_atac_z)
cell_order = cell_order[hc_res2$order]
##cell_order = sort(cell_order)
cell_order = factor(cell_order, levels = annot_order)
cell_order_df = data.frame(cell_order)
colnames(cell_order_df) = "annotation"

## Colors
##hm_colors1 = rev(hcl.colors(n=5, palette=col_pal))[col_range]
##hm_colors1 = colorRampPalette(colors = hm_colors1)(100)
hm_colors1 = rev(hcl.colors(100, "RdBu"))
##hm_colors2 = rev(hcl.colors(n=5, palette=col_pal))[col_range]
##hm_colors2 = colorRampPalette(colors = hm_colors2)(100)
hm_colors2 = rev(hcl.colors(100, "Spectral"))

## Plot
##pdf(file=paste0(out_dir, out_name, ".pdf"), width = 8, height = 12)
png(filename = paste0(out_dir, out_name, ".png"), width = 800, height = 1200)
pheatmap(pb_atac_z[hc_res1$order, hc_res2$order],
         color = hm_colors1,
         cluster_rows = F, cluster_cols = F,
         labels_col = "", labels_row = "",
         annotation_col = cell_order_df, 
         annotation_colors = list(annotation = color_scheme),
         annotation_names_row = F, 
         annotation_names_col = F)
dev.off()

##pdf(file=paste0(out_dir, out_name2), width = 8, height = 12)
png(filename = paste0(out_dir, out_name2, ".png"), width = 800, height = 1200)
pheatmap(pb_rna_z[hc_res1$order, hc_res2$order], 
         color = hm_colors2,
         cluster_rows = F, cluster_cols = F, 
         show_rownames = F, show_colnames = F,
         annotation_col = cell_order_df, 
         annotation_colors = list(annotation = color_scheme),
         annotation_names_row = F, 
         annotation_names_col = F)
dev.off()

#################
### Figure 4d ###
#################
df = data.frame(CRE_genes)
colnames(df) = c("gene", "count")
ggplot(df, aes(count)) + 
  geom_histogram(binwidth=1, color="black", fill="deepskyblue") +
  xlim(0, 50) +
  labs(x="Number of linked CREs", y="Frequency") +
  geom_vline(xintercept = 8, color="red", size=1, linetype="dashed") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20))

#################
### Figure 4e ###
#################
df = cor_res_expr_peak %>% filter(r > 0.2 & pvalue_adj < 0.05)
nCREs = table(df$gene)
nCREs = data.frame(nCREs)
colnames(nCREs) = c("gene", "nLinks")
nCREs$gene = as.character(nCREs$gene)
nCREs = nCREs[order(nCREs$nLinks), ]
nCREs$idx = 1:nrow(nCREs)
N = 22
cond = nCREs$idx > (nrow(nCREs) - N)
nCREs$label= ifelse(cond, nCREs$gene, "")
nCREs$highlight = ifelse(cond, "yes", "no")

p = ggplot(nCREs, aes(x=idx, y=nLinks, color=highlight, label=label)) + 
  geom_point(size=1.5) + 
  geom_text_repel(size=4) + 
  scale_color_manual(values = c("yes"="red", "no"="lightgrey")) +
  labs(x="Rank", y="Number of linked CREs") +
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=20),
        legend.position = "None")
p + ylim(c(0, 50))

#################
### Figure 4f ###
#################
a = names(nCREs_sub)
b = egenes
ggvenn(
  list("GPC" = a, "SE-linked gene" = b),
  show_percentage = F,
  stroke_size = 1, set_name_size = 6, text_size = 8
) 

#################
### Figure 4g ###
#################
axis_breaks = c(-1, 0, 1)
axis_limits = c(-1.5, 1.5)

ggplot(df, aes(Expr, Activity, color=Score, label=Top)) + 
  geom_point() + 
  geom_text_repel(size=3, color="red") +
  scale_color_gradient2(low="blue", mid="white", high="red") +
  scale_x_continuous(breaks = axis_breaks, limits = axis_limits) +
  scale_y_continuous(breaks = axis_breaks, limits = axis_limits) +
  labs(x="Expression (Z-score)", y="Regulon activity (Z-score)", 
       color="Chromvar score") +
  facet_wrap(.~Subpopulation, nrow = 1, scales="free") + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size=20),
        axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))

#################
### Figure 4h ###
#################
gene = "SEZ6L"

## Identify CREs linked to target gene
peaks = rownames(top_CREs[top_CREs$gene == gene, ])
peaks_gr = peaksToGranges(peaks)
chr = unlist(strsplit(peaks[1], split="-"))[1]

## Identify tss 
gene_tss = tss.hg19[tss.hg19$gene == gene, "tss"]

## Identify target gene and its CREs in the closest link df 
gene_peak = closest_features[closest_features$gene_name == gene & closest_features$peaks %in% peaks, ]
##peaks_gr = peaks_gr[peaks_gr$peaks %in% gene_peak$peaks, ]

## Add start, end, mid point of the link and associated TSS 
gene_peak$start = start(peaks_gr)
gene_peak$end = end(peaks_gr)
gene_peak$tss = sapply(rownames(gene_peak),
                       function(x){
                         res = NULL
                         min_dist = 1E9
                         for (i in gene_tss){
                           tmp = min(abs(gene_peak[x, "start"] - i), 
                                     abs(gene_peak[x, "end"] - i))
                           if (tmp < min_dist){
                             res = i
                             min_dist = tmp
                           }
                         }
                         return(res)
                       })
gene_peak$mid = floor((gene_peak$start+gene_peak$end)/2)

## Generate a Grange object for mid point of each CRE and its associated TSS and add that as a Link obj to seurat obj
peak_link = sapply(rownames(gene_peak), 
                   function(x){
                     a = min(gene_peak[x, "tss"], gene_peak[x, "mid"])
                     b = max(gene_peak[x, "tss"], gene_peak[x, "mid"])
                     paste(chr, a, b, sep="-")
                   })
peak_link_gr = peaksToGranges(peak_link)
peak_link_gr$score = top_CREs[gene_peak$peaks, "r"]
Links(seurat_obj) <- peak_link_gr

## overlap between peaks and SEs
overlap_se = findOverlaps(peaks_gr, all_ses) 
all_ses[all_ses$peak_idx %in% subjectHits(overlap_se), ]

## Gene identified by integrating multiomics
gene_tss_w_CREs = as.numeric(names(table(gene_peak$tss)))
chr_region = paste(chr, (min(gene_tss_w_CREs)-1E5), (max(gene_tss_w_CREs)+1E5), sep="-")

## Coverage plot
Idents(seurat_obj) = seurat_obj$annotation_sep
cov_plot = CoveragePlot(
  object = seurat_obj,
  region = gene,
  #extend.upstream = 1E5,
  #extend.downstream = 1E5,
  annotation = F, 
  peaks = F, 
  links = F
)
cov_plot = cov_plot & scale_fill_manual(values=color_scheme)

## gene plot
gene_plot <- AnnotationPlot(
  object = seurat_obj,
  region = chr_region
)

## peak coordinates 
peak_plot = PeakPlot(
  object = seurat_obj,
  region = chr_region,
  peaks = peaks_gr
)

## Links
link_plot <- LinkPlot(
  object = seurat_obj,
  region = chr_region
) 
link_plot = link_plot & scale_color_gradientn(colors=brewer.pal(9, "Reds")) & labs(color="correlation\ncoefficient")

CombineTracks(
  plotlist = list(cov_plot, peak_plot, gene_plot, link_plot),
  heights = c(10, 1, 2, 3)
)
