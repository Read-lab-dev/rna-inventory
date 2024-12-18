library(slingshot)
library(destiny)
library(Biobase)
rm(list = ls())
gsc.seu <-qs::qread("gsc.seu.qs")
gsc.seu <- subset(gsc.seu,celltype=="tCYC",invert=T)
gsc.seu <- subset(gsc.seu,orig.ident=="Eng_DMSO")
dat <- as.data.frame(t(GetAssayData(gsc.seu,layer = "count"))) %>% tibble::rownames_to_column(var = "Cell")
ct <- as.ExpressionSet(dat)
ct$type <- as.vector(gsc.seu$celltype)
##Diffusion map
# housekeepers <- c('ACTB', 'GAPDH')  
# normalizations <- colMeans(exprs(ct)[housekeepers,])
# ct_norm <- ct
# exprs(ct_norm)<- expr(ct_norm)-normalizations
dm <- DiffusionMap(ct,n_pcs = 50)
palette(cube_helix(6)) # configure color palette
plot(dm, 1:2,legend_main="type")
qs::qsave(dm,"dm.qs")
library(ggplot2)
DC <- eigenvectors(dm)[,1:2];colnames(DC) <- c("DC_1","DC_2");rownames(DC) <- colnames(gsc.seu)
gsc.seu[["DC"]] <- CreateDimReducObject(embeddings = DC, key = "DC_")
DimPlot(gsc.seu,reduction = "DC")

#######SlingShot#########
gsc.sce <- as.SingleCellExperiment(gsc.seu)
gsc.sce <- slingshot(gsc.sce, reducedDim = 'DC',clusterLabels="celltype")
# start.clus="tOPC",end.clus=c("tAC","tNPC") # no clusters

pdf(file = "./slingshot/time_DMSO.pdf",height = 4,width = 5)
colors <- rev(met.brewer('Hiroshige',50))
plot(reducedDims(gsc.sce)$DC, col = colors[cut(gsc.sce$slingPseudotime_1,breaks=50)],
     pch=16,axes=FALSE)
lines(SlingshotDataSet(gsc.sce), lwd=2,col="black")
legend.col(colors,seq(0,0.5,length.out=50))
dev.off()

pdf(file = "./slingshot/cluster_DMSO.pdf",height = 4,width = 5)
plot(reducedDims(gsc.sce)$DC, col = brewer.pal(4,'Set1')[gsc.sce$celltype],
     pch=16,axes=FALSE)
lines(SlingshotDataSet(gsc.sce), lwd=2,col="black")
legend('right', legend = levels(gsc.sce$celltype)[-2], col = brewer.pal(4,'Set1')[-2], pch = 16, bty = 'n')
dev.off()

# slingshot_df <- data.frame(colData(gsc.sce)[, names(colData(gsc.sce)) != 'slingshot', drop=FALSE])
# # Plot Slingshot pseudotime vs cell stage. 
# ggplot(slingshot_df, aes(x = slingPseudotime_1, y = celltype, 
#                                              colour =celltype))+
#   ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
#   scale_color_tableau() + theme_classic() +
#   xlab("Slingshot pseudotime") + ylab("Timepoint") +
#   ggtitle("Cells ordered by Slingshot pseudotime")

legend.col <- function(col, lev) {
  opar <- par()
  n <- length(col)
  bx <- par("usr")
  box.cx <- c(bx[1] + (bx[2] - bx[1]) / 1000, bx[1] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[4])
  box.sy <- (bx[4] - bx[3]) / n
  xx <- rep(box.cx, each = 2)
  par(xpd = TRUE)
  for (i in 1:n) {
    yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + (box.sy * i), box.cy[1] + (box.sy * i), box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
  }
  par(new = TRUE)
  plot(0, 0, type = "n", ylim = c(min(lev), max(lev)), yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
  axis(side = 2, las = 2, tick = FALSE, line = 0.25)
  par <- opar
}
