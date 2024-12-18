# Ensure Giotto Suite and GiottoData packages are installed
if(!"Giotto" %in% installed.packages()) {
  devtools::install_github("drieslab/Giotto@suite")
}

if(!"Giotto" %in% installed.packages()) {
  devtools::install_github("drieslab/GiottoData")
}

library(Giotto)
library(GiottoData)

# Ensure the Python environment for Giotto has been installed
genv_exists = checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to install the Giotto environment.
  installGiottoEnvironment()
}

# to automatically save figures in save_dir set save_plot to TRUE
temp_dir = getwd()
myinstructions = createGiottoInstructions(save_dir = temp_dir,
                                          save_plot = TRUE,
                                          show_plot = TRUE)

# Provide path to visium folder
data_path = paste0('./DMG_visium/data/DMG1_spaceranger_out')

# Create Giotto object
visium_lungcancer = createGiottoVisiumObject(visium_dir = data_path,
                                             expr_data = 'filter',
                                             png_name = 'tissue_lowres_image.png',
                                             gene_column_index = 2,
                                             instructions = myinstructions)

# check metadata
pDataDT(visium_lungcancer)

# check available image names
showGiottoImageNames(visium_lungcancer) # "image" is the default name

# show aligned image
spatPlot(gobject = visium_lungcancer, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7)

visium_lungcancer <- filterGiotto(gobject = visium_lungcancer,
                                  expression_threshold = 1,
                                  feat_det_in_min_cells = 50,
                                  min_det_feats_per_cell = 1000,
                                  expression_values = 'raw',
                                  verbose = T)
visium_lungcancer <- normalizeGiotto(gobject = visium_lungcancer, scalefactor = 6000, verbose = T)
visium_lungcancer <- addStatistics(gobject = visium_lungcancer)

spatPlot2D(gobject = visium_lungcancer, show_image = T, point_alpha = 0.7)

spatPlot2D(gobject = visium_lungcancer, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_feats', color_as_factor = F)
visium_lungcancer <- calculateHVF(gobject = visium_lungcancer)
visium_lungcancer <- runPCA(gobject = visium_lungcancer)
screePlot(visium_lungcancer, ncp = 30)
plotPCA(gobject = visium_lungcancer)
visium_lungcancer <- runUMAP(visium_lungcancer, dimensions_to_use = 1:10)
plotUMAP(gobject = visium_lungcancer)

# Create shared nearest network (SNN) and perform leiden clustering
visium_lungcancer <- createNearestNetwork(gobject = visium_lungcancer, dimensions_to_use = 1:10, k = 30)
visium_lungcancer <- doLeidenCluster(gobject = visium_lungcancer, spat_unit = 'cell', feat_type = 'rna', resolution = 0.4, n_iterations = 1000)

# visualize UMAP cluster results
plotUMAP(gobject = visium_lungcancer, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2)

# visualize expression and spatial results
spatPlot2D(gobject = visium_lungcancer, cell_color = 'leiden_clus')

spatDimPlot(gobject = visium_lungcancer, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 2, dim_show_legend = T, spat_show_legend = T, spat_point_size = 2)

# Cell type marker detection
# Gini markers
gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_lungcancer,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_featss = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)

# get top 2 genes per cluster and visualize with violin plot
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats
violinPlot(visium_lungcancer, feats = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right')

# cluster heatmap
plotMetaDataHeatmap(visium_lungcancer,
                    selected_feats = topgenes_gini,
                    metadata_cols = c('leiden_clus'),
                    x_text_size = 10, y_text_size = 10)

# Cell type marker detection
# Scran markers
scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_lungcancer,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')

# get top 2 genes per cluster and visualize with violin plot
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats
violinPlot(visium_lungcancer, feats = unique(topgenes_scran),
           cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right')

#########PAGE###############
Giotto_SC <- createGiottoObject(expression = "STARmap_3D_data_expression.txt",
                                spatial_locs ="STARmap_3D_data_cell_locations.txt")

# umap plots
# Create PAGE matrix
# PAGE matrix should be a binary matrix with each row represent a gene marker and each column represent a cell type
# There are several ways to create PAGE matrix
# 1.1 create binary matrix of cell signature genes
# small example #
Tcells_markers = c("CD2", "CD3D", "CD3E", "CD3G")
macrophage_markers = c("MARCO", "CSF1R", "CD68", "GLDN", "APOE", "CCL3L1", "TREM2", "C1QB", "NUPR1", "FOLR2", "RNASE1", "C1QA")
dendritic_markers = c("CD1E", "CD1C", "FCER1A", "PKIB", "CYP2S1", "NDRG2")
mast_markers = c("CMA1", "TPSAB1", "TPSB2")
Bcell_markers = c("IGLL5", "MZB1", "JCHAIN", "DERL3", "SDC1", "MS$A1", "BANK1", "PAX5", "CD79A")
Bcell_PB_markers = c("PRDM1", "XSP1", "IRF4")
Bcell_mem_markers = c("MS4A1", "IRF8")
housekeeping_markers = c("ACTB", "GAPDH", "MALAT1")
neutrophils_markers = c("FCGR3B", "ALPL", "CXCR1", "CXCR2", "ADGRG3", "CMTM2", "PROK2", "MME", "MMP25", "TNFRSF10C")
pdcs_markers = c("SLC32A1", "SHD", "LRRC26", "PACSIN1", "LILRA4", "CLEC4C", "DNASE1L3", "SCT", "LAMP5")

signature_matrix = makeSignMatrixPAGE(sign_names = c('T_Cells', 'Macrophage', 'Dendritic', 'Mast', 'B_cell', 'Bcell_PB', 'Bcells_memory',
                                                     'Housekeeping', 'Neutrophils', 'pDCs'),
                                      sign_list = list(Tcells_markers,
                                                       macrophage_markers,
                                                       dendritic_markers,
                                                       mast_markers,
                                                       Bcell_markers,
                                                       Bcell_PB_markers,
                                                       Bcell_mem_markers,
                                                       housekeeping_markers,
                                                       neutrophils_markers,
                                                       pdcs_markers))

# 1.3 enrichment test with PAGE

markers_scran = findMarkers_one_vs_all(gobject=Giotto_SC, method="scran",
                                       expression_values="normalized", cluster_column = "Class", min_feats=3)

top_markers <- markers_scran[, head(.SD, 10), by="cluster"]
celltypes<-levels(factor(markers_scran$cluster))
sign_list<-list()
for (i in 1:length(celltypes)){
  sign_list[[i]]<-top_markers[which(top_markers$cluster == celltypes[i]),]$feats
}

PAGE_matrix_3 = makeSignMatrixPAGE(sign_names = celltypes,
                                   sign_list = sign_list)

#  runSpatialEnrich() can also be used as a wrapper for all currently provided enrichment options
visium_lungcancer = runPAGEEnrich(gobject = visium_lungcancer, sign_matrix = signature_matrix, min_overlap_genes = 1)

# 1.4 heatmap of enrichment versus annotation (e.g. clustering result)
cell_types = colnames(signature_matrix)
plotMetaDataCellsHeatmap(gobject = visium_lungcancer,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types,
                         spat_enr_names = 'PAGE',
                         x_text_size = 8,
                         y_text_size = 8,
                         show_plot = T,
                         save_param = list(save_name="7_a_metaheatmap"))