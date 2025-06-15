setwd("~/rna/shinyR")
# devtools::install_github('SGDDNB/ShinyCell')
library(ShinyCell)
library(shiny)
seu = qs::qread("../BT109VP/gsc.seu.qs")
seu$cycling <- ifelse(seu$cycscore1>0,"cycling","non-cycling")
save(cycling)
scConf = createConfig(seu)
# Delete excessive metadata and rename some metadata
scConf = delMeta(scConf,colnames(seu@meta.data)[c(6:10,12,14:17)])
scConf = modMetaName(scConf, 
                     meta.to.mod = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt",
                                     "log10GenesPerUMI","celltype","H3K27M"), 
                     new.name = c("Treatment", "nUMI","nFeature","Mito.Percent",
                                  "log10GenesPerUMI","Cell.State","H3K27M-detection"))
# Modify colours and labels
seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
scConf = modColours(scConf, meta.to.mod = "Treatment", 
                    new.colours=  c("#CF5A79","#D77B5A","#5DA373"))
scConf = modColours(scConf, meta.to.mod = "Cell.State", 
                    new.colours=  c("#DF9ED4","#544799","#5DA373","#D77B5A"))
scConf = modColours(scConf, meta.to.mod = "cycling", 
                    new.colours=  c("darkred","#8080804D"))
scConf = modColours(scConf, meta.to.mod = "H3K27M-detection", 
                    new.colours=  c("darkblue","#8080804D"))
scConf = modLabels(scConf, meta.to.mod = "Treatment", 
                   new.labels = c("DMSO", "LongVP", "ShortVP"))

showLegend(scConf)
scConf = reorderMeta(scConf, scConf$ID[c(1,6,8,2,3,5,4,7)])
showOrder(scConf)
scConf = modDefault(scConf, "Cell.State", "Treatment")

makeShinyApp(seu, scConf, gene.mapping = TRUE,
             shiny.title = "BT109 Engrafted Tumor Cells",
             shiny.footnotes = "Created by Zhengang Hu",
             default.gene1 = "EGFR", default.gene2 = "YAP1",
             default.dimred = c("dim21","dim22"))

runApp("bt109/", host ="0.0.0.0",launch.browser=T ,port = 1111)


##########RSconnect#########

#https://www.shinyapps.io/

rsconnect::setAccountInfo(name='read',
                          token='xxx',
                          secret='xxx')
rsconnect::deployApp('bt109/')

rsconnect::setAccountInfo(name='read', token='xxx',
                          secret='xxx')

rsconnect::deployApp('bt109/')
