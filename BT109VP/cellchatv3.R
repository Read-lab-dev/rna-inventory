library(Seurat)
library(CellChat)
options('future.globals.maxSize'=1024*1024^3)
# gsc.seu <- DietSeurat(gsc.seu)
# aa <- subset(int.seu,celltype=="eGSC"&gfp=="eGFP-")
gsc.seu <- qs::qread("gsc.seu.qs")
int.seu <- qs::qread("int.seu.NEW1.qs")
gsc.seu$celltype <- paste0(gsc.seu$celltype,"-like")
Idents(gsc.seu) <- gsc.seu$celltype
aa.seu <- subset(int.seu,celltype=="eGSC",invert=T)
aa.seu <- merge(gsc.seu,aa.seu)
Idents(aa.seu) <- aa.seu$celltype

DMSO.seu <- subset(aa.seu,orig.ident=="Eng_VP_L")

cellchat <- createCellChat(object = DMSO.seu, group.by = "ident", assay = "RNA")
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
vertex.receiver = c(1:4)# a numeric vector.
vertex.sender = c(5:6)
cellchat@idents
df.net <- subsetCommunication(cellchat,targets.use = 1:3,sources.use =4:13)
qs::qsave(cellchat,file = 'cellchat-new-L.qs')
