rm(list = ls())
gc()
###########SC-RNAseq############
## We load the required packages
library(Seurat)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
gsc.seu <-qs::qread("../gsc.seu.qs")
net <- get(data("dorothea_hs", package = "dorothea"))
net
# Extract the normalized log-transformed counts
mat <- as.matrix(gsc.seu@assays$RNA@data)
# Run wmean
acts <- run_wmean(mat=mat, net=net, .source='tf', .target='target',
                  .mor='mor', times = 100, minsize = 5)
acts

# Extract norm_wmean and store it in tfswmean in pbmc
gsc.seu[['tfswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = gsc.seu) <- "tfswmean"
# Scale the gsc.seu
gsc.seu <- ScaleData(gsc.seu)
gsc.seu@assays$tfswmean@data <- gsc.seu@assays$tfswmean@scale.data

DefaultAssay(object = gsc.seu) <- "tfswmean"
p1 <- DimPlot(gsc.seu, reduction = "dim2", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(gsc.seu,features = c("YAP1"),reduction = "dim2") & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('YAP1 activity')
p3 <- FeaturePlot(gsc.seu, features = c("YAP1"),reduction = "dim2") + ggtitle('YAP1 expression')
DefaultAssay(object = gsc.seu) <- "tfswmean"
p1 | p2 | p3

n_tfs <- 25
# Extract activities from object as a long dataframe
df <- t(as.matrix(gsc.seu@assays$tfswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(gsc.seu)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long gsc.seu frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

# Plot
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
