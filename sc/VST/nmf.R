rm(list = ls())
setwd("~/rna/sc/VST")
mye.seu <- qs::qread("mye_new.seu.qs")
library(NMF)
DefaultAssay(mye.seu) <- "RNA"
mye.seu = NormalizeData(mye.seu,verbose=F)
mye.seu = FindVariableFeatures(mye.seu,selection.method = "vst",nfeatures= 2000,verbose=F)
mye.seu = ScaleData(mye.seu,verbose=F)
# mat <- as.matrix(mye.seu@assays$RNA@counts[rownames(mye.seu@assays$RNA@scale.data),])
mat <- mye.seu@assays$RNA@scale.data
mat[mat<=0] <- 0
mat <- mat[rowSums(mat)>0,]
res.rank <- NMF::nmf(mat,
                     rank = 6:10,
                     .opt="vP10",
                     method = "snmf/r")
plot(res.rank)

nmf.rank4 <- NMF::nmf(mat,
                      rank = 4,
                      nrun=50,
                      .opt="vP10",
                      method = "snmf/r")
fs <- extractFeatures(nmf.rank4, 50L)
fs <- lapply(fs, function(x) rownames(nmf.rank4)[x])
fs <- do.call("rbind",fs)
rownames(fs) <- paste0("cluster", 1:4)
DT::datatable(t(fs))

index <- extractFeatures(nmf.rank4,50L) 
sig.order <- unlist(index)
mat <- mye.seu@assays$RNA@scale.data
NMF.Exp.rank4 <- mat[sig.order,]
NMF.Exp.rank4 <- na.omit(NMF.Exp.rank4) #sig.order有时候会有缺失值
group <- predict(nmf.rank4) # 提出亚型
NMF.Exp.rank4 <- NMF.Exp.rank4[,order(group)]
# NMF.Exp.rank4 <- scale(NMF.Exp.rank4,center = T,scale = T)
NMF.Exp.rank4[NMF.Exp.rank4>3] <-3
NMF.Exp.rank4[NMF.Exp.rank4<-3] <--3
ac=data.frame(NMF=sort(group))
ac$NMF <- paste0("C",ac$NMF)
jco <- list(NMF=c(C1="#2874C5",C2="#EABF00",C3="#C6524A",C4="limegreen"))
pheatmap::pheatmap(NMF.Exp.rank4,cluster_cols=F,cluster_rows = F,
                   show_rownames = F,show_colnames = F,
                   annotation_col = ac,
                   annotation_colors = jco,height = 4,width = 5,
                   filename = "./NMF/Heatmap_NMF3.pdf")
lineage <- as.data.frame(t(fs))
DoHeatmap(mye.seu,features =rownames(NMF.Exp.rank4),group.by = "cluster",label = F)
ggsave("./NMF/heatmap.pdf")
save(nmf.rank4,file = "./NMF/nmf.rank4.Rds")
colnames(lineage) <- c("N","AC","OPC","cyc")
write.csv(lineage,file="./NMF/module4.csv")
#########
s.f =1:4
cell1 <- colnames(mye.seu)
cell2 <- colnames(coef(nmf.rank4))
cells <- intersect(cell1, cell2)
mye.seu <- mye.seu[,cells]
mye.seu <- RunPCA(mye.seu, verbose = F)
mye.seu@reductions$nmf <- mye.seu@reductions$pca
mye.seu@reductions$nmf@cell.embeddings <- t(coef(nmf.rank4)[,cells])    
colnames(mye.seu@reductions$nmf@cell.embeddings) <- c("PC_1","PC_2","PC_3","PC_4")
mye.seu@reductions$nmf@feature.loadings <- basis(nmf.rank4)  
Embeddings(mye.seu,reduction = "nmf")
mye.seu <- RunTSNE(mye.seu,reduction='nmf', dims=1:4)
mye.seu <- RunUMAP(mye.seu,reduction='nmf', dims=1:4)
mye.seu$cluster <- predict(nmf.rank4)
Idents(mye.seu) <- mye.seu$cluster
new.cluster.ids <- c("Neu-like","AC-like","OPC","Cycling")
names(new.cluster.ids) <- levels(mye.seu)
mye.seu <- RenameIdents(mye.seu, new.cluster.ids)
DimPlot(mye.seu,reduction = "tsne")
qs::qsave(mye.seu,"mye.seu.qs")
#################

library(RcppML)
# cells is a Seurat object
model <- RcppML::nmf(mye.seu@assays$RNA@data, k = 30, verbose = F, seed = 1234)
W <- model$w # Components Matrix
H <- model$h # Activity Matrix

# Name genes/components for latent factor matrix
colnames(W) <- paste0("component", 1:30)
rownames(W) <- rownames(mye.seu) # Gene names

# Name components/cells for activity matrix
colnames(H) <- colnames(mye.seu) # Cell IDs
rownames(H) <- paste0("component", 1:30)
# This code assumes the above code block has been already run
require(dplyr) #used for %>% piping
N <- 50
top_genes <- list()
for (i in 1:dim(W)[2]) {
  name <- colnames(W)[i]
  top_genes[[name]] <- sort(W[, i], decreasing = T) %>% head(N) %>% names()
}
