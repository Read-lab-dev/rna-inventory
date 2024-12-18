setwd("~/rna2/cri")
library(Seurat)
int.seu <- qs::qread("gsc.seu.qs")
lineage <- read.csv("inflammatory",header = F,na.strings = "")
int.seu <- AddModuleScore(int.seu,features = list(na.omit(lineage$V1)),nbin = 30,ctrl = 100,name = "inf")

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis,color=inf))+
  geom_point(alpha=1,size=0.5)+
  scale_colour_distiller(palette="RdBu")+
geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))
FeaturePlot(int.seu,features = "inf1")+scale_colour_distiller(palette="RdBu")

library(AUCell) 
library(clusterProfiler) 
inflam <- list(read.table("inflammatory")$V1)
names(inflam)<- "inflammatory"
cells_rankings <- AUCell_buildRankings(int.seu@assays$RNA@data, nCores=10, plotStats=TRUE)

c5 <- read.gmt("/home/hzg/rna/MsigDB/c5.go.bp.v2023.1.Hs.symbols.gmt")
genename <- c("GOBP_ASTROCYTE_DIFFERENTIATION","GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS",
              "GOBP_OLIGODENDROCYTE_DEVELOPMENT","GOBP_MESENCHYMAL_STEM_CELL_DIFFERENTIATION",
              "GOBP_NEURON_DEVELOPMENT","GOBP_ACUTE_INFLAMMATORY_RESPONSE","GOBP_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE"
              )
genename <- c("GOBP_MESENCHYMAL_STEM_CELL_PROLIFERATION","GOBP_POSITIVE_REGULATION_OF_MESENCHYMAL_STEM_CELL_PROLIFERATION",
              "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION","GOBP_MESENCHYMAL_CELL_PROLIFERATION","GOBP_MESENCHYMAL_CELL_DIFFERENTIATION")

genename <- "GOBP_CELL_CYCLE_CHECKPOINT_SIGNALING"
geneSets <- lapply(genename, function(x){c5$gene[c5$term == x]})
names(geneSets) <- genename
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

set.seed(123)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
DimPlot(int.seu,cells.highlight = cells_assignment$GOBP_POSITIVE_REGULATION_OF_MESENCHYMAL_STEM_CELL_PROLIFERATION$assignment,sizes.highlight = 0.1)
int.seu$GOBP_MESENCHYMAL_CELL_DIFFERENTIATION <- as.numeric(getAUC(cells_AUC)["GOBP_MESENCHYMAL_CELL_DIFFERENTIATION", ])
int.seu$GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION <- as.numeric(getAUC(cells_AUC)["GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION", ])
VlnPlot(int.seu,features = "GOBP_MESENCHYMAL_CELL_DIFFERENTIATION",split.by = "orig.ident",group.by = "state")
FeaturePlot(int.seu,features = c("GOBP_MESENCHYMAL_CELL_DIFFERENTIATION","GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"))&viridis::scale_color_viridis(option="C")&theme(title = element_text(size = 10))

int.seu$HALLMARK_P53_PATHWAY <- as.numeric(getAUC(cells_AUC)["HALLMARK_INFLAMMATORY_RESPONSE",])

FeaturePlot(int.seu,features = c("HALLMARK_P53_PATHWAY"))&viridis::scale_color_viridis(option="F")&theme(title = element_text(size = 10))

VlnPlot(int.seu,features = "TP53",group.by = "orig.ident")

auc_data <- cbind(int.seu@meta.data[,1:26],t(cells_AUC@assays@data$AUC))

library(ggpubr)
ggboxplot(auc_data, x = "orig.ident", y = "cycling.score", fill = "orig.ident",
          color = "black", palette = c("#CF5A79","#D2C564","#5DA373"),
          xlab = "Treatment")+stat_compare_means(comparisons =list(c("CrizL","DMSO"),c("CrizS", "DMSO")))+facet_wrap(vars(state))

c2 <- read.gmt("/home/hzg/rna/MsigDB/c2.all.v2023.1.Hs.symbols.gmt")
genename <- c("VERHAAK_GLIOBLASTOMA_PRONEURAL","VERHAAK_GLIOBLASTOMA_NEURAL",
              "VERHAAK_GLIOBLASTOMA_CLASSICAL","VERHAAK_GLIOBLASTOMA_MESENCHYMAL")
geneSets <- lapply(genename, function(x){c2$gene[c2$term == x]})
names(geneSets) <- genename
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc_data <- cbind(auc_data[,1:26],t(cells_AUC@assays@data$AUC))

set.seed(10242023)
ggboxplot(auc_data, x = "orig.ident", y = "VERHAAK_GLIOBLASTOMA_CLASSICAL", fill = "orig.ident",
          color = "black", palette = c("#CF5A79","#D2C564","#5DA373"),
          xlab = "Treatment")+stat_compare_means(comparisons =list(c("CrizL","DMSO"),c("CrizS", "DMSO")))+facet_wrap(vars(state))

cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE,nCores = 10) 
DimPlot(int.seu,cells.highlight = cells_assignment$VERHAAK_GLIOBLASTOMA_CLASSICAL$assignment,sizes.highlight = 0.1)
