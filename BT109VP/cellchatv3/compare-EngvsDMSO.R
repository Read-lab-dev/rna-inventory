rm(list = ls())
gc()
setwd("~/rna/BT109VP/cellchatv3")
library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
cellchat.DMSO <- qs::qread("cellchat-new-DMSO.qs")
# cellchat.DMSO <- subsetCellChat(cellchat.DMSO,idents.use = c("tNPC","tAC","tOPC","tCYC"),invert = T)
cellchat.org <- qs::qread("cellchat-new-Org.qs")
object.list <- list(Org = cellchat.org, Engrafted = cellchat.DMSO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
#interaction.num.strength.VPs.png

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# (D) Circle plot showing the differential number of interactions or interaction strength among coarse cell types
group.cellType <- factor(c(rep("Tumor", 4),"Astroglia","Astroglia","Cycling",rep("Neuron",4 )))
object.list$Org <- mergeInteractions(object.list$Org, group.merged =as.character(group.cellType[-c(1:4)]))
object.list$Engrafted <- mergeInteractions(object.list$Engrafted, group.merged =group.cellType)

object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#Compare the major sources and targets in a 2D space

#(A) Identify cell populations with significant changes in sending or receiving signals

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(unlist(num.link)), max(unlist(num.link)))

object.list <- lapply(object.list, function(x) {netAnalysis_computeCentrality(x)})

gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
