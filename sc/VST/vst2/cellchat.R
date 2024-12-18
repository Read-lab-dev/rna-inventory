rm(list = ls())
gc()
library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
ptm = Sys.time()
options(stringsAsFactors = FALSE)
options(future.globals.maxSize=2000*1024^2)
output_dir <-"./cellchatv2"
int.seu <- qs::qread("sch.seu.qs")
dir.create(output_dir)
setwd(output_dir)
Idents(int.seu) <- int.seu$celltype
cellchat <- createCellChat(object = int.seu, group.by = "ident", assay = "RNA")
cellchatDB <- CellChatDB.human 
showDatabaseCategory(cellchatDB)
# Show the structure of the database
dplyr::glimpse(cellchatDB$interaction)

CellChatDB.use <- cellchatDB
# use all cellchatDB for cell-cell communication analysis
# cellchatDB.use <- cellchatDB # simply use the default cellchatDB. We do not suggest to use it in this way because cellchatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling).
# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multicore", workers = 10) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)

ptm = Sys.time()
cellchat <- computeCommunProb(cellchat,type = "triMean")

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
qs::qsave(cellchat,file = "cellchat.qs")


cellchat <- qs::qread("cellchat.qs")
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = c(7)# a numeric vector.

vertex.sender = c(2:6,8:12)
df.net <- subsetCommunication(cellchat,targets.use = vertex.receiver)
pathways.show <- c("ADGRG","GAS","MK","NOTCH") 
par(mfrow=c(2,4),xpd=TRUE)
for (i in 1:length(pathways.show)) {
  netVisual_aggregate(cellchat, signaling = pathways.show[i])
}
for (i in 1:length(pathways.show)) {
  netVisual_aggregate(cellchat, signaling = pathways.show[i],targets.use = vertex.receiver)
}
# Circle plot
pathways.show <- c("MK") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
netAnalysis_contribution(cellchat, signaling = pathways.show)

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "circle",out.format = "svg")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.png"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

netVisual_chord_gene(cellchat, targets.use = 11, legend.pos.x = 15)

netAnalysis_contribution(cellchat.DMSO, signaling = "NOTCH")

plotGeneExpression(cellchat, signaling = "NOTCH", split.by = "datasets", colors.ggplot = T, type = "violin")

VlnPlot(int.seu,features = "NOTCH1",split.by = "orig.ident")


########NPC signaling####
cellchat <- qs::qread("cellchat-VP-S.qs")
signaling <- c("NRXN","NRG","NCAM","CDH","MK","NTS","PTN","CDH")
netVisual_chord_gene(cellchat, sources.use = 13, lab.cex = 0.4,signaling = signaling,
                     title.name = "NPC-like cells")
signaling <- c("MK","PTN","NRXN","PTPR","CNTN")
netVisual_chord_gene(cellchat, targets.use = 13, signaling = signaling,lab.cex = 0.4)
