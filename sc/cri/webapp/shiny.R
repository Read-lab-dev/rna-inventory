setwd("/home/hzg/rna/sc/cri/web/")
# devtools::install_github('SGDDNB/ShinyCell')
library(ShinyCell)
library(shiny)
seu = qs::qread("./DMSO-tumor.qs")
seu 
# select metadata and rename some metadata
seu@meta.data <- seu@meta.data[,c(1:4,7,20,23,25,26)]
scConf = createConfig(seu)
scConf = modMetaName(scConf, 
                     meta.to.mod = scConf$ID, 
                     new.name = c("Ident", "nUMI","nFeature","Mito.Percent",
                                  "cluster","cycling","Cell.State","ZM-fusion","MET-ex14"))
# Modify colours and labels
seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
scConf = modColours(scConf, meta.to.mod = "Ident", 
                    new.colours=  c("#D77B5A"))
scConf = modColours(scConf, meta.to.mod = "Cell.State", 
                    new.colours=  c("#544799","#CF5A79","#D2C564","#5DA373"))
scConf = modColours(scConf, meta.to.mod = "cycling", 
                    new.colours=  c("darkred","#8080804D"))
scConf = modColours(scConf, meta.to.mod = "ZM-fusion", 
                    new.colours=  c("darkblue","#8080804D"))
scConf = modColours(scConf, meta.to.mod = "MET-ex14", 
                    new.colours=  c("#924099","#8080804D"))

showLegend(scConf)
scConf = reorderMeta(scConf, scConf$ID[c(7,6,8,9,2,3,4,5,1)])
showOrder(scConf)
scConf = modDefault(scConf, "Cell.State", "cluster")

makeShinyApp(seu, scConf, gene.mapping = TRUE,
             shiny.title = "BT85 Engrafted Tumor Cells",
             shiny.footnotes = "Created by Zhengang Hu",
             default.gene1 = "MET", default.gene2 = "PTPRZ1",
             default.dimred = c("Dim21","Dim22"))

runApp("BT85/", host ="0.0.0.0",launch.browser=T ,port = 1111)

setwd("/home/hzg/rna/sc/cri/webapp/")
##########RSconnect#########
rsconnect::setAccountInfo(name='read', 
                          token='7DC568EBA96673A27C55D2FB1A401328', 
                          secret='Cwhag8n26Gy6mdBJFReu0izyPXNl6zalB+6opyMg')
rsconnect::deployApp('BT85/')
