
int.seu[["tf"]]<-CreateAssayObject(counts=getAUC(sub_regulonAUC))
DefaultAssay(int.seu) <- "tf"
int.seu <- FindVariableFeatures(int.seu)
int.seu <- NormalizeData(int.seu)
int.seu <- ScaleData(int.seu)
int.seu <- RunPCA(int.seu)
int.seu <- RunUMAP(int.seu,dims = 1:30)
int.seu <-FindNeighbors(int.seu)

DimPlot(int.seu, group.by = "celltype",reduction = "umap", label = T,repel = T,label.size = 4,label.box = T,cols = MetBrewer::met.brewer("Klimt",12),label.color = "grey100",
        pt.size = 0.1)+
  theme_void()+
  theme(legend.position = "none",title = element_blank())

##############
library(Seurat)
library(tidyverse)
library(magrittr)
visuNetwork <- function(regulon.name,n){
  adj.ls <- regulon.name %>% map(~{
    tmp <- adj[which(adj$TF == regulon.name), ]
    tmp <- tmp[order(tmp$importance, decreasing = T),]
    loci <- unique(c(1:n, grep(regulon.name, tmp$target)))
    tmp <- tmp[loci,]
    return(tmp)
  })
  adj.sub <- Reduce('rbind', adj.ls)
  ## generate network
  edge.df <- adj.sub
  colnames(edge.df) <- c('from', 'to', 'weights')
  edge.df$from <- as.character(edge.df$from)
  edge.df$to <- as.character(edge.df$to)
  vertex <- unique(c(edge.df$from, edge.df$to))
  net <- graph_from_data_frame(d = edge.df, vertices = vertex, directed = T)
  vsize <- scale(edge.df$weights, center = min(quantile(edge.df$weights)), scale = max(quantile(edge.df$weights)) - min(quantile(edge.df$weights)))
  plot(
    net,
    # vertex.label = vlabel,
    vertex.label.cex = 1,
    vertex.label.dist = 1,
    vertex.label.color = 'black',
    vertex.size = c(1, vsize+0.1)*10,
    #vertex.size = vsize*10,
    vertex.frame.color = NA,
    edge.arrow.size = 0.5,
    #  vertex.color = vcl, 
    main = regulon.name
  )
}
adj <- read.table('scenic/all data/adj.sample.tsv',header = T)
names(regulons) <- gsub("[(+)]", "", names(regulons))
lapply("TEAD4", n=20,visuNetwork)
ggsave("TEAD4.PDF")
