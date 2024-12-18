rm(list = ls())
setwd("~/rna/BT109VP/VIPER")
library(Seurat)
BiocManager::install("bcellViper")
BiocManager::install("viper")
devtools::install_github('saezlab/dorothea')
library(viper)
library(dorothea)
gsc.seu <- qs::qread("../gsc.seu.qs")
load("~/rna/BT109VP/VIPER/gbm_u133a_cleaner_mas5_tfregulon.rda")
## We read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

gsc.seu <- run_viper(gsc.seu, regulon,
                 options = list(method = "scale", minsize = 4, 
                                eset.filter = FALSE, cores = 4,  
                                verbose = FALSE))

DefaultAssay(object = gsc.seu) <- "dorothea"
gsc.seu.markers <- FindAllMarkers(object = gsc.seu, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(gsc.seu.markers)
