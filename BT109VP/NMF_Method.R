###NMF##
# devtools::install_github("zdebruine/RcppML")
# devtools::install_github("zdebruine/singlet")
library(singlet)
library(Seurat)
library(dplyr)

lineage <- as.data.frame(t(fs))
dat <- as.data.frame(nmf.rank4@fit@W)
dat$symbol <- rownames(dat)
{ library(clusterProfiler)
  library(dplyr)
  
  ##Creating Gene List for GSEA
  dat2 <- dat[,c(4,5)]
  colnames(dat2)[1] <- "log2FoldChange"
  genelist <- bitr(dat2$symbol,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(dat2, by= c("SYMBOL"="symbol")) %>% mutate(rank = rank(log2FoldChange,ties.method = "random")) %>%
    arrange(desc(rank))
  
  gsea_input <- genelist$log2FoldChange
  
  names(gsea_input) <- genelist$ENTREZID
  
  ##Select your Gene Set
  
  dir='~/backup/Bulk_Analysis/MsigDB/'
  
  gmts <- list.files(dir,pattern = 'gmt')
  
  gmts
  
  #Start GSEA Analysis
  library(GSEABase)
  
  gsea_results <- lapply(gmts, function(gmtfile){
    # gmtfile=gmts[2]
    geneset <- read.gmt(file.path(dir,gmtfile)) 
    print(paste0('Now process the ',gmtfile))
    egmt <- GSEA(gsea_input, TERM2GENE=geneset, verbose=T,nPermSimple = 100000)
    head(egmt)
    return(egmt)
  })
  gsea_results[[1]] <- setReadable(gsea_results[[1]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  
  gsea_results[[2]] <- setReadable(gsea_results[[2]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  
  gsea_results_df <- do.call(rbind, gsea_results_list)
}
library(enrichplot)
enrichplot::gseaplot2(gsea_results[[2]],c("GOBP_NEURON_FATE_COMMITMENT","GOBP_AXON_DEVELOPMENT","GOBP_NEURON_DEVELOPMENT"),
                      color = seurat_colors,pvalue_table = T)

save(gsea_results,gsea_results_df,file = "GSEA_Eng_VP_cyc.Rdata")
write.csv(gsea_results_df,file = "GSEA_Eng_VP_cyc.csv")
qs::qsave(gsc.seu,file = "gsc.nmf.qs")


gsc.seu <- qs::qread("gsc.seu.qs")

library(NMF)
DefaultAssay(gsc.seu) <- "RNA"
gsc.seu = NormalizeData(gsc.seu,verbose=F)
gsc.seu = FindVariableFeatures(gsc.seu,selection.method = "vst",nfeatures= 3000,verbose=F)
gsc.seu = ScaleData(gsc.seu,verbose=F)
mat <- gsc.seu@assays$RNA@scale.data
mat <- as.matrix(gsc.seu@assays$RNA@scale.data)
mat[mat<=0] <- 0
mat <- mat[rowSums(mat)>0,]
res.rank <- NMF::nmf(mat,
                     rank = 2:10,
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
mat <- gsc.seu@assays$RNA@scale.data
NMF.Exp.rank4 <- mat[sig.order,]
NMF.Exp.rank4 <- na.omit(NMF.Exp.rank4) #sig.order有时候会有缺失值
group <- predict(nmf.rank4) # 提出亚型
NMF.Exp.rank4 <- NMF.Exp.rank4[,order(group)]
# NMF.Exp.rank4 <- scale(NMF.Exp.rank4,center = T,scale = T)
NMF.Exp.rank4[NMF.Exp.rank4>3] <-3
NMF.Exp.rank4[NMF.Exp.rank4<-0] <--3
ac=data.frame(NMF=sort(group))
ac$NMF <- paste0("C",ac$NMF)
jco <- list(NMF=c(C1="#2874C5",C2="#EABF00",C3="#C6524A",C4="limegreen"))
pheatmap::pheatmap(NMF.Exp.rank4,cluster_cols=F,cluster_rows = F,
                   show_rownames = F,show_colnames = F,
                   annotation_col = ac,
                   annotation_colors = jco,height = 4,width = 5,
                   filename = "./Heatmap_NMF4.pdf")
lineage <- as.data.frame(t(fs))
DoHeatmap(gsc.seu,features =rownames(NMF.Exp.rank4),group.by = "cluster",label = F)
ggsave("./NMF/heatmap.pdf")
save(nmf.rank4,file = "./NMF/nmf.rank4.Rds")
colnames(lineage) <- c("N","AC","OPC","cyc")
write.csv(lineage,file="./NMF/module4.csv")
#########
s.f =1:4
cell1 <- colnames(gsc.seu)
cell2 <- colnames(coef(nmf.rank4))
cells <- intersect(cell1, cell2)
gsc.seu <- gsc.seu[,cells]
gsc.seu <- RunPCA(gsc.seu, verbose = F)
gsc.seu@reductions$nmf <- gsc.seu@reductions$pca
gsc.seu@reductions$nmf@cell.embeddings <- t(coef(nmf.rank4)[,cells])    
colnames(gsc.seu@reductions$nmf@cell.embeddings) <- c("PC_1","PC_2","PC_3","PC_4")
gsc.seu@reductions$nmf@feature.loadings <- basis(nmf.rank4)  
Embeddings(gsc.seu,reduction = "nmf")
gsc.seu <- RunTSNE(gsc.seu,reduction='nmf', dims=1:4)
gsc.seu <- RunUMAP(gsc.seu,reduction='nmf', dims=1:4)
gsc.seu$cluster <- predict(nmf.rank4)
Idents(gsc.seu) <- gsc.seu$cluster
new.cluster.ids <- c("Neu-like","AC-like","OPC","Cycling")
names(new.cluster.ids) <- levels(gsc.seu)
gsc.seu <- RenameIdents(gsc.seu, new.cluster.ids)
DimPlot(gsc.seu,reduction = "tsne")
qs::qsave(gsc.seu,"gsc.seu.qs")
#################
### Figure 2i ###
#################
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gghighlight)
library(MetBrewer)
library(RColorBrewer)
library(tibble)
gsc.seu <- qs::qread("gsc.seu.qs")
ns.seu <- qs::qread("neurosphere.qs")

seurat_colors <- as.character(met.brewer("Klimt", 8))
colors_continuous <- as.character(rev(met.brewer("Hiroshige", 100)))
colors_scale <- as.character(rev(met.brewer("Johnson",100)))
sub_cols <- c("#CF5A79","#544799")

lineage <- read.csv("./NMF/module4.csv",header = T,na.strings = "")

gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$NPC)),nbin = 30,ctrl = 100,name = "NPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$AC)),nbin = 30,ctrl = 100,name = "AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$OPC)),nbin = 30,ctrl = 100,name = "OPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$cyc)),nbin = 30,ctrl = 100,name = "cycscore")
set.seed(20230829)

lineage.score <- pmax(gsc.seu$NPC1,gsc.seu$AC1)
lineage.score.plot <- lineage.score
lineage.class <- ifelse(gsc.seu$NPC1>gsc.seu$AC1,"NPC","AC")
# lineage.score.plot[lineage.score.plot<0] <- runif(length(lineage.score.plot[lineage.score.plot<0]),0,0.15)
lineage.data <- data.frame(lineage.score,lineage.score.plot,lineage.class)
lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"] <- -lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"]
lineage.data$lineage.class <- ifelse(gsc.seu$OPC1>pmax(gsc.seu$AC1,gsc.seu$NPC1),"OPC",lineage.class)
# lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"]),-0.1,0.1)
# 
# lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"]),-0.1,0.1)
lineage.data$stemness <- gsc.seu$OPC1-lineage.data$lineage.score
lineage.data <- data.frame(lineage.data,gsc.seu@meta.data)
lineage.data$cycling <- ifelse(lineage.data$cycscore1>0,"cycling","non-cycling")

lineage.data$eng <- ifelse(lineage.data$orig.ident%in%c("BT109_DMSO","BT109_VP"),"unengrafted","engrafted")
lineage.data$vp <- ifelse(lineage.data$orig.ident%in%c("Eng_VP_L","Eng_VP_S","BT109_VP"),"VP_Treated","DMSO")
lineage.data$orig.ident <- factor(lineage.data$orig.ident,levels = c("Eng_DMSO","Eng_VP_S","Eng_VP_L"))
seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
sub_cols <- c("#CF5A79","#544799")
lineage.data$SOX2 <- GetAssayData(gsc.seu)["SOX2",]

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=orig.ident))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors[c(1,3,5,7,9)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))
ggsave("./2dim-VP.pdf",height = 3,width = 7.5)

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=celltype))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values =seurat_colors[c(3,7,2,5)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "top",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )+
ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=cycling))+
  geom_point(size=0.8)+scale_colour_manual(values = c("darkred","#8080804D"))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "top",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )+
  ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=H3K27M))+
  geom_point(size=0.8)+scale_colour_manual(values = c("darkblue","#8080804D"))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "top",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )

ggsave("./2dim-class.pdf",height = 4,width = 8)

lineage.data$SOX2 <- GetAssayData(gsc.seu)["SOX2",]
ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=SOX2))+
  geom_point(alpha=1,size=0.8)+scale_color_gradientn(colors=c("grey" ,"red"))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )
ggsave("./2dim-sox2.pdf",height = 4,width = 6.5)

qs::qsave(gsc.seu,file = "gsc.seu.qs")
write.csv(lineage.data,file = "lineage.data.csv")
#############Neurosphere#########
gsc.seu <- qs::qread("neurosphere.qs")
DefaultAssay(gsc.seu) <-"RNA"
gsc.seu <- merge(gsc.seu,qs::qread("gsc.seu.qs"))

lineage <- read.csv("./NMF/module4.csv",header = T,na.strings = "")

gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$NPC)),nbin = 30,ctrl = 100,name = "NPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$AC)),nbin = 30,ctrl = 100,name = "AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$OPC)),nbin = 30,ctrl = 100,name = "OPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$cyc)),nbin = 30,ctrl = 100,name = "cycscore")
set.seed(20230829)

lineage.score <- pmax(gsc.seu$NPC1,gsc.seu$AC1)
lineage.score.plot <- lineage.score
lineage.class <- ifelse(gsc.seu$NPC1>gsc.seu$AC1,"NPC","AC")
# lineage.score.plot[lineage.score.plot<0] <- runif(length(lineage.score.plot[lineage.score.plot<0]),0,0.15)
lineage.data <- data.frame(lineage.score,lineage.score.plot,lineage.class)
lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"] <- -lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"]
lineage.data$lineage.class <- ifelse(gsc.seu$OPC1>pmax(gsc.seu$AC1,gsc.seu$NPC1),"OPC",lineage.class)
# lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="AC"&lineage.data$lineage.score>0,"lineage.score.plot"]),-0.1,0.1)
# 
# lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"] <-runif(length(lineage.data[lineage.data$lineage.class=="NPC"&lineage.data$lineage.score<0,"lineage.score.plot"]),-0.1,0.1)
lineage.data$stemness <- gsc.seu$OPC1-lineage.data$lineage.score
lineage.data <- data.frame(lineage.data,gsc.seu@meta.data)
lineage.data$cycling <- ifelse(lineage.data$cycscore1>0,"cycling","non-cycling")

lineage.data$eng <- ifelse(lineage.data$orig.ident%in%c("BT109_DMSO","BT109_VP"),"unengrafted","engrafted")
lineage.data$vp <- ifelse(lineage.data$orig.ident%in%c("Eng_VP_L","Eng_VP_S","BT109_VP"),"VP_Treated","DMSO")
lineage.data$orig.ident <- factor(lineage.data$orig.ident,levels = c("Eng_DMSO","Eng_VP_S","Eng_VP_L"))
seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
sub_cols <- c("#CF5A79","#544799")

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=orig.ident))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors[c(1,3,5,7,9)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))
ggsave("./2dim-VP-ns.pdf",height = 3,width = 7.5)

ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=cycling))+
  geom_point(size=0.8)+scale_colour_manual(values = c("darkred","#8080804D"))+
  theme_few()+xlab("AC-like <--------------> NPC-like")+
  theme(
    legend.position = "top",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )+
  ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=lineage.class))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors[c(1,3,5,7,9)])+
  theme_few()+xlab("AC-like <--------------> NPC-like")+
  theme(
    legend.position = "top",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  )

ggsave("./2dim-class-ns.pdf",height = 4,width = 6.5)


###############Jessa##################
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$NPC)),nbin = 30,ctrl = 100,name = "NPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$AC)),nbin = 30,ctrl = 100,name = "AC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$OPC)),nbin = 30,ctrl = 100,name = "OPC")
gsc.seu <- AddModuleScore(gsc.seu,features = list(na.omit(lineage$cyc)),nbin = 30,ctrl = 100,name = "cycscore")
set.seed(20230829)

lineage.score <- pmax(gsc.seu$NPC1,gsc.seu$AC1)
lineage.score.plot <- lineage.score
lineage.class <- ifelse(gsc.seu$NPC1>gsc.seu$AC1,"NPC","AC")
lineage.data <- data.frame(lineage.score,lineage.score.plot,lineage.class)
lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"] <- -lineage.data[lineage.data$lineage.class=="AC","lineage.score.plot"]
lineage.data$lineage.class <- ifelse(gsc.seu$OPC1>pmax(gsc.seu$AC1,gsc.seu$NPC1),"OPC",lineage.class)
lineage.data$stemness <- gsc.seu$OPC1-lineage.data$lineage.score
lineage.data <- data.frame(lineage.data,gsc.seu@meta.data)
lineage.data$cycling <- ifelse(lineage.data$cycscore1>0,"cycling","non-cycling")
seurat_colors <- met.brewer("Klimt",13)
lineage.data <- lineage.data[lineage.data$cell!="P-6253_S-8498_TGTACAGAGTTTCTTC-1",]
ggplot(data = lineage.data,aes(x=lineage.score.plot,y=stemness,color=orig.ident))+
  geom_point(alpha=1,size=0.8)+scale_colour_manual(values = seurat_colors)+
  theme_few()+xlab("AC-like <--------------> NPC-like")+gghighlight(unhighlighted_params =list(colour="grey90",alpha=0.8))+
  theme(
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid")
  ) +facet_wrap(vars(orig.ident))

gsc.seu$deg <- paste0(gsc.seu$orig.ident,"_",gsc.seu$celltype)

Idents(gsc.seu) <- gsc.seu$celltype
aa= FindAllMarkers(gsc.seu,test.use = "MAST",logfc.threshold = 1)
##########Add reduction########
dim2 <- lineage.data[,c(2,4)]
identical(rownames(dim2),colnames(gsc.seu))
colnames(dim2) <- c("dim2_1","dim2_2")
dim2$dim2_1 <- as.numeric(dim2$dim2_1)
dim2$dim2_2 <- as.numeric(dim2$dim2_2)
dim2_emb<- as.matrix(dim2)
gsc.seu[["dim2"]] <- CreateDimReducObject(dim2_emb,key = "dim2_")
DimPlot(gsc.seu,reduction = "dim2")

######COUNT########
count_dat <- lineage.data %>% group_by(orig.ident,celltype) %>% summarise(Count=n()) %>% mutate(Frac = Count / sum(Count) * 100) %>% as.data.frame()

count_dat$orig.ident <- factor(count_dat$orig.ident,levels = c("Eng_DMSO","Eng_VP_S","Eng_VP_L"))
ggplot(data=count_dat, mapping=aes(x=celltype,y=Frac,fill=orig.ident))+
  geom_bar(stat="identity",width=0.5,position='dodge')+
  scale_fill_manual(values=c('#999999','#E69F00',"#544799","#CF5A79"))+
  geom_text(stat='identity',aes(y=Frac+2,label=paste0(round(Frac,0),"%")),color="black",position = position_dodge(0.5),fontface="bold")+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

######Density######
library(ggthemes)
library(ggplot2)
library(ggnewscale)
hierarchy.result <- gsc.seu@meta.data
hierarchy.result <- cbind(hierarchy.result,gsc.seu@reductions$dim2@cell.embeddings)
hierarchy.result$YAP1 <- GetAssayData(gsc.seu)["YAP1",]
hierarchy.result$TEAD1 <- GetAssayData(gsc.seu)["TEAD1",]
ggplot(data = hierarchy.result,aes(x=dim2_1,y=dim2_2))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  theme_few()+xlab("AC-like <--------------> NPC-like")+ylab("stemness")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = hierarchy.result[hierarchy.result$YAP1>0,],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+ggtitle("BT109 Engrafted Tumor Cells")+scale_fill_viridis_c(option = "B",name="YAP1+Cells")+scale_alpha_continuous(range = c(0,1),guide=guide_none())

ggplot(data = hierarchy.result,aes(x=dim2_1,y=dim2_2,color=SOX2))+
  geom_point(alpha=1,size=0.5)+
  theme_few()+xlab("AC-like <--------------> NPC-like")+ylab("stemness")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+scale_color_viridis_c(option = "D")

lineage.data$YAP1 <-GetAssayData(gsc.seu)["EGFR",]

library(ggridges)
ggplot(lineage.data,aes(x=YAP1,y=celltype,fill=stat(x)))+
  geom_density_ridges_gradient( jittered_points = TRUE,
                                position = position_points_jitter(width = 0.05, height = 0),
                                point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7)+  scale_fill_viridis_c(option = "E",name="")+
  labs(title = 'stemness')+
  theme_ridges()

#################Heatmap##################
dat <- nmf.rank4@consensus

pheatmap::pheatmap(dat,cluster_cols=F,cluster_rows = F,
                   show_rownames = F,show_colnames = F,
                   annotation_col = ac,
                   annotation_colors = jco,height = 4,width = 5,
                   filename = "./NMF/Heatmap_NMF10.pdf")
dev.new()
par(mfrow=c(1,2))
consensusmap(nmf.rank4)
library(Seurat)
library(ggplot2)
seurat_colors <- c("#DF9ED4","#CF5A79","#D77B5A","#D2C564","#5DA373","#406E89","#544799", "#924099")
gsc.seu <- qs::qread("gsc.seu.qs")
fs <- extractFeatures(nmf.rank4,50L)
fs <- lapply(fs, function(x) rownames(nmf.rank4)[x])
fs <- sapply(fs,"[",i=1:max(sapply(fs,length)))
fs <- fs[,c(2,1,3,4)]
DoHeatmap(gsc.seu,features = fs,draw.lines	=F)+ 
  scale_fill_gradientn(colors = c("#64a5c5", "white", "#b4212a"))

dat <- gsc.seu@assays$RNA$scale.data[fs,]
dat <- cor(t(dat))

row_annot = data.frame(
  meta = rep(c("NPC", "AC", "OPC","CYC"), times = c(50, 50, 50, 50)),
  row.names = rownames(dat))
ann_colors = list(meta=c(NPC="#D77B5A", AC="#544799", OPC="#D2C564",CYC="#CF5A79"))
breaksList = seq(-1,1, by = 0.1)

pheatmap::pheatmap(dat,cluster_cols=F,cluster_rows = F,
                   annotation_row =row_annot,annotation_col = row_annot,
                   annotation_colors = ann_colors,
                   breaks = breaksList,gaps_col = c(50,100,150),gaps_row = c(50,100,150),
                   color = colorRampPalette(c("#64a5c5", "white", "#b4212a"))(length(breaksList)),
                   show_rownames = F,show_colnames = F,
                   filename = "heatmap-metaprogram.pdf",height = 6,width = 7)


