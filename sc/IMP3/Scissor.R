#Scissor
# devtools::install_github("sunduanchen/Scissor")
library(scAB)
library(Scissor)
library(Seurat)
library(ggthemes)
library(ggplot2)
library(ggnewscale)
setwd("~/rna/sc/IMP3")
source("./bulk/scissor_code.R")
bulk_dataset <- expr_norm
  
phenotype <- ifelse(metadata$Treat=="sh75",1,0)
tag <- c('ctrl', 'shIMP3')
  
sc_dataset <- run_seurat(gbm.seu,verbose = FALSE)
infos4 <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag,alpha = 0.05,
                  cutoff = 0.2,family = "binomial", Save_file = "Scissor_GBM301_sh75_alpha0.051.RData",Load_file = "Scissor_GBM301_sh75_alpha0.05.RData")

infos4 <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag,alpha = 0.1,
                  cutoff = 0.2,family = "binomial", Save_file = "Scissor_GBM301_sh3_alpha0.11.RData",Load_file = "Scissor_GBM301_sh3_alpha0.1.RData")

Scissor_select <- rep("unselcted", ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos4$Scissor_pos] <- "shRIOK2"
Scissor_select[infos4$Scissor_neg] <- "ctrl"
gbm.seu <- AddMetaData(gbm.seu, metadata = Scissor_select, col.name = "scissor")
# gbm.seu$cellstate <- hierarchy.result$cellstate
UMAP_celltype <- DimPlot(gbm.seu, reduction ="dim2",label = T,pt.size=0.0001,group.by = "cellstate")+scale_color_manual(values = c("#CF5A79","#544799","#D2C564","#5DA373"))+theme_few()&NoLegend()

UMAP_scissor <- DimPlot(gbm.seu, reduction = 'dim2', group.by = 'scissor', cols = c('royalblue','darkred','grey'), pt.size=0.0001)+theme_few()
UMAP_celltype+UMAP_scissor
ggsave(UMAP_celltype+UMAP_scissor,file="Scissor_GBM301_shRIOK2_select.pdf",height = 3.5,width = 8)

DimPlot(gbm.seu, reduction = 'dim2', group.by = 'scissor', cols = c('royalblue','darkred','grey'), pt.size=0.0001,split.by = "Sample",ncol = 5)+theme_few()+
  theme(strip.background = element_rect(
    color="black", fill="#FC4E07", size=1.5, linetype="solid"
  )
)

hierarchy.result$scissor <-Scissor_select
p1=ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$scissor=="ctrl",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+ggtitle("GBM301_ctrl")+scale_fill_viridis_c()+scale_alpha_continuous(range = c(0,1),guide=guide_none())
p2=ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$scissor=="shRIOK2",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+ggtitle("GBM301_shRIOK2")+scale_fill_viridis_c()+scale_alpha_continuous(range = c(0,1),guide=guide_none())
p1+p2
ggsave(p1+p2,file="Scissor_GBM301_shRIOK2.pdf",height = 5,width = 12)

ggplot(data = hierarchy.result,aes(x=X,y=Y))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$scissor=="shRIOK2",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+ggtitle("GBM301 on Neftel et.al")+scale_fill_viridis_c(option = "B",name="shRIOK2")+scale_alpha_continuous(range = c(0,1),guide=guide_none())+
  new_scale_fill()+
  new_scale("alpha")+
  stat_density2d(data = hierarchy.result[hierarchy.result$scissor=="ctrl",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+scale_fill_viridis_c(option = "D",name="ctrl")+scale_alpha_continuous(range = c(0,1),guide=guide_none())

ggsave(file="Scissor_GBM301_shRIOK2_ctrl.pdf",height = 4,width = 6)

numbers <- length(infos4$Scissor_pos) + length(infos4$Scissor_neg)
result1 <- reliability.test(X, Y, network, alpha = 0.1, family = "binomial", cell_num = numbers,n=10,nfold = 2)

Idents(gbm.seu) <- gbm.seu$scissor

aa <-FindMarkers(gbm.seu,ident.1 = "shRIOK2",ident.2 = "ctrl",test.use = "MAST",logfc.threshold = 0)

resdata <- aa
resdata$gene <- rownames(resdata)
colnames(resdata)[2] <- "log2FoldChange"
write.csv(resdata,file = "DEseq2_sh75_scissor.csv")

