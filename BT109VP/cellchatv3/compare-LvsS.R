rm(list = ls())
gc()
library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
cellchat.DMSO <- qs::qread("cellchat-new-S.qs")
cellchat.DMSO <- subsetCellChat(cellchat.DMSO,idents.use = c("tNPC-like","tAC-like"),invert = T)
cellchat.S <- qs::qread("cellchat-new-L.qs")
object.list <- list(VP.S = cellchat.DMSO, VP.L = cellchat.S)
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

group.cellType <- c(rep("eGSC", 2), rep("Neuron",4 ),"AC","AC","Cycling")
# group.cellType <- factor(group.cellType, levels = c("AC", "N", "Tumor"))
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

#(B) Identify the signaling changes of specific cell populations
#Furthermore, we can identify the specific signaling changes of Inflam.DC and cDC1 between NL and LS.
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "tOPC-like",comparison = c(1, 2))
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "tCYC-like",comparison = c(1, 2))
patchwork::wrap_plots(plots = list(gg1,gg4))
ggsave(wrap_plots(plots = list(gg1,gg4)),file="VP.L.change.pdf",height = 5,width = 12)
#Identify signaling groups based on their functional similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional",comparison=c(1,2))
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional",comparison=c(1,2))
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional",comparison=c(1,2))
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5,comparison=c(1,2))

rankSimilarity(cellchat, type = "functional")

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE,comparison = c(1,2))
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE,comparison = c(1,2))

gg1 + gg2

ggsave(gg1+gg2,file="inforflow-VP.S.png")

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 20, height = 24)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 20, height = 24)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 20, height = 24, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 20, height = 24, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 20, height = 24, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 20, height = 24, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#"overall.signaling.patterns.VP.s.png"
##Identify dysfunctional signaling by comparing the communication probabities

netVisual_bubble(cellchat,sources.use = 7, targets.use = c(1:4),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:4),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in VP.S", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:4),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in VP.S", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

#Identify dysfunctional signaling by using differential expression analysis
pos.dataset = "VP.S"
features.name = paste0(pos.dataset, ".merged")

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "VP.S",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "VP.S",ligand.logFC = -0.05, receptor.logFC = NULL)

table(net.down$pathway_name)
df <- findEnrichedSignaling(object.list[[2]], features = c("NOTCH1"), idents = c("T-OPC"), pattern ="outgoing")

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 1, targets.use = c(11:14), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 1, targets.use = c(11:14), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

# visualize the enriched ligands in the first condition
par(mfrow = c(1,2), xpd=TRUE)
computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)

pathways.show <- c("NRXN") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show,targets.use = 13, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CNTN") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("NRXN") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("DMSO", "VP.S")) # set factor level
plotGeneExpression(cellchat, signaling = "BMP", split.by = "datasets", colors.ggplot = T, type = "violin")

save(object.list, file = "cellchat_object.list_VP.S.RData")
save(cellchat, file = "cellchat_merged_VP.S.RData")
