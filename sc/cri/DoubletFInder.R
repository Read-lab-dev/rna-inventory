## DoubletFinder
library(DoubletFinder)
library(Seurat)
library(dplyr)

sweep.res.list <- paramSweep_v3(int.seu, PCs = 1:20, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = 0.039                     # 5000细胞对应的doublets rate是3.9%
homotypic.prop <- modelHomotypic(int.seu$seurat_clusters)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(int.seu)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
int.seu <- doubletFinder_v3(int.seu, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, 
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)


DimPlot(int.seu, reduction = "umap", group.by = "DF.classifications_0.25_0.2_2009")

int.seu<- subset(int.seu,subset = DF.classifications_0.25_0.2_2009 =="Singlet")
int.seu <- standard10X(int.seu, nPCs=20, res=0.8)
DimPlot(int.seu,cells.highlight = WhichCells(int.seu,expression = eGFP==0))

org.seu <- subset(int.seu,eGFP==0)
org.seu <- standard10X(org.seu, nPCs=20, res=0.5)

gsc.seu <- subset(int.seu,eGFP==0&seurat_clusters=="8",invert = TRUE)
gsc.seu <- standard10X(gsc.seu, nPCs=20, res=0.5)

DimPlot(org.seu,label = T,group.by = "RNA_snn_res.0.6")+DimPlot(gsc.seu,label = T,group.by = "RNA_snn_res.0.6")

qs::qsave(org.seu,file = "org.seu.qs")
qs::qsave(gsc.seu,file = "gsc.seu.qs")

DimPlot(org.seu,group.by = "RNA_snn_res.0.6",label = T,label.box = T,cols = cluster.col)
DimPlot(gsc.seu,group.by = "RNA_snn_res.0.6",label = T,label.box = T,cols = cluster.col)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

gsc.seu <- CellCycleScoring(gsc.seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
org.seu <- CellCycleScoring(org.seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(gsc.seu,group.by = "Phase",label = T,label.box = T,cols = sample(MetBrewer::met.brewer("Klimt",6),3))
DimPlot(org.seu,group.by = "Phase",label = T,label.box = T,cols = sample(MetBrewer::met.brewer("Klimt",6),3))

Idents(org.seu)<- org.seu$RNA_snn_res.0.6
features <- list(
  Astrocyte = c("HES1","GFAP","AGT","GLUL"),
  Neuroepithelial = c("KLF4","VIM","VCAN"),
  Neuron    = c("RBFOX3","DCX","MAP2"),
  pNPC      = c("PAX6","SOX2","GLI3","TOP2A"),
  EarlyNeuroectodermProgenitors  = c("PCP4","NEFL","NEFM"),
  Fib       = c("COL1A1","PDGFRA","DCN"),
  GliomaStemCell= c("EGFR","CD44","PROM1","FUT4","NEU1","NES","eGFP","MET"))
DotPlot(int.seu,features = features)

all.markers <- FindAllMarkers(org.seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25) 

top10.markers <- all.markers[-c(grep(c("^RP[SL]"),rownames(all.markers)),
                                grep(c("^MT-"),rownames(all.markers))),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

cellmarker <- NULL
for (i in levels(org.seu)) {
  tmp <-top10.markers %>% filter(cluster==i) %>% dplyr::select(gene) %>% as.data.frame()
  
  cluster <- paste(tmp$gene,collapse = ",")
  
  cellmarker <- rbind(cellmarker,cluster)
}
rownames(cellmarker) <- paste0("cluster",levels(org.seu))
write.csv(cellmarker,file = "cellmarker_org.csv")

cluster.col <- sample(MetBrewer::met.brewer("Klimt",20),20)

cell.prop<-as.data.frame(prop.table(table(gsc.seu$RNA_snn_res.0.6, gsc.seu$orig.ident)))
colnames(cell.prop)<- c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  scale_fill_manual(values = cluster.col)+
  geom_bar(stat="identity",position="fill")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))+coord_flip()

int.seu$state <- gsc.module$state
cell.prop<-as.data.frame(prop.table(table(int.seu$state, int.seu$orig.ident)))
colnames(cell.prop)<- c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  scale_fill_manual(values = cluster.col)+
  geom_bar(stat="identity",position="fill")+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))+coord_flip()


