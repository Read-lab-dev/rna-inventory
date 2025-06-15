setwd("/home/hzg/rna/sc/tarun")
# devtools::install_github('SGDDNB/ShinyCell')
library(ShinyCell)
library(shiny)
seu = qs::qread("./int.seu.qs")
scConf = createConfig(seu)
# Delete excessive metadata and rename some metadata
scConf = delMeta(scConf,colnames(seu@meta.data)[c(6:8,10,11)])
scConf = modMetaName(scConf, 
                     meta.to.mod = scConf$ID, 
                     new.name = c("Cell.Line", "nUMI","nFeature","Mito.Percent",
                                  "log10GenesPerUMI","Cell.Cluster","Cell.Phase","Cell.State","Cell.State.Hybrid"))
# Modify colours and labels
seurat_colors <- c("#3D5495","#DF9ED4","#CC4E56","#924099","#DF9D5D","#D2698D","#624699","#6AA771","#C5C166","#428181")
scConf = modColours(scConf, meta.to.mod = "Cell.Line", 
                    new.colours=  c("#CF5A79","#3D5495","#5DA373","#924099"))
scConf = modColours(scConf, meta.to.mod = "Cell.State", 
                    new.colours=  c("#544799","#CF5A79","#D2C564","#5DA373"))
scConf = modColours(scConf, meta.to.mod = "Cell.State.Hybrid", 
                    new.colours=  seurat_colors[1:7])
scConf = modColours(scConf, meta.to.mod = "Cell.Cluster", 
                    new.colours=  seurat_colors)

showLegend(scConf)
scConf = reorderMeta(scConf, scConf$ID[c(1,6,8,9,2,3,5,4,7)])
showOrder(scConf)
scConf = modDefault(scConf, "Cell.State", "Cell.Line")

setwd("./webapp/")
makeShinyApp(seu, scConf, gene.mapping = TRUE,
             shiny.title = "GBM Neurospheres Lines",
             shiny.footnotes = "Created by Zhengang Hu",
             default.gene1 = "EGFR", default.gene2 = "YAP1",
             default.dimred = c("dim21","dim22"),
             shiny.dir = "gbmnsp/")
runApp("gbmnsp/", host ="0.0.0.0",launch.browser=T ,port = 1111)


##########RSconnect#########
rsconnect::setAccountInfo(name='read',
                          token='xxx',
                          secret='xxx')
rsconnect::deployApp('gbmnsp/')
