rm(list = ls())
gsc.seu <- readRDS("~/rna/scRNA_analysis/BT109_VP/gsc.seu.Rds")
state.color <- c("#544799","#CF5A79","#D2C564","#5DA373")
######Prop#########

gsc.seu <- Seurat::CellCycleScoring(gsc.seu,s.features = cc.genes.updated.2019$s.genes,
                         g2m.features =cc.genes.updated.2019$g2m.genes )
cell.prop <- as.data.frame(prop.table(table(gsc.seu$orig.cycling, gsc.seu$state)))
colnames(cell.prop) <- c("cluster", "sample", "proportion")
ggplot(cell.prop, aes(sample, proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = state.color)+
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.1, "cm"),
    legend.key.size = unit(10, "pt"),
    axis.text = element_text(size = 5, face = "bold"),
    title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_legend(title = NULL))
ggsave("cell.prop_hybrid.pdf", width = 4, height = 3)

cell.prop <- as.data.frame(prop.table(table(gsc.seu$state, gsc.seu$orig.ident)))

DimPlot(gsc.seu,group.by = "state")

colnames(cell.prop) <- c("cluster", "sample", "proportion")

ggplot(cell.prop, aes(sample, proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = state.color)+
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.1, "cm"),
    legend.key.size = unit(10, "pt"),
    axis.text = element_text(size = 5, face = "bold"),
    title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_legend(title = NULL))

#########
gsc.module <- data.frame(gsc.seu$MES11,gsc.seu$MES21,gsc.seu$AC1,
                         gsc.seu$OPC1,gsc.seu$NPC11,gsc.seu$NPC21)
colnames(gsc.module) <- c("MES1","MES2","AC","OPC","NPC1","NPC2")
gsc.module$state <- apply(gsc.module,1,which.max)
gsc.module$state <- gsub(1,"MES1",gsc.module$state)
gsc.module$state <- gsub(2,"MES2",gsc.module$state)
gsc.module$state <- gsub(3,"AC",gsc.module$state)
gsc.module$state <- gsub(4,"OPC",gsc.module$state)
gsc.module$state <- gsub(5,"NPC1",gsc.module$state)
gsc.module$state <- gsub(6,"NPC2",gsc.module$state)
gsc.module$state.value <- apply(gsc.module[,1:6],1,max)
gsc.module$state.hybrid <- gsc.seu$state.hybrid

gsc.module <- gsc.module %>% group_by(state) %>% arrange(state.hybrid,state,desc(state.value))

ggplot(data = gsc.module,aes(x=xaxis,y=yaxis))+
  geom_point(alpha=1,size=0.5,color="grey90")+
  geom_vline(xintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  geom_hline(yintercept=0, linetype="dashed", colour="grey", linewidth=0.5) + 
  theme_few()+xlab("AC-like <--------------> MES-like")+ylab("AC-like <--------------> OPC-like")+
  theme(axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "gray60", size = 1, linetype = "solid"))+
  stat_density2d(data = gsc.module[gsc.module$state.hybrid=="non",],
                 aes(fill=..density..,alpha=..density..),geom ="raster",contour = F)+
  scale_fill_gradientn(colors = colors_continuous)


heatmap.mat <- as.matrix(gsc.module[,c(3,1,2,5,6,4)])
n <- t(scale(heatmap.mat))
n[n>3] <- 3
n[n<-3] <- -3
colgroup <- as.data.frame(gsc.module$state.hybrid)
colnames(colgroup) <- "celltype"
rownames(colgroup) <- paste0("V",1:ncol(n))
colnames(n) <- paste0("V",1:ncol(n))
anno_col_color <- list(
  celltype=c(cycling="#CF5A79",hybrid="#544799",non="#d2c564")
)
break_point <- c(table(gsc.module$state.hybrid)[1],
                 table(gsc.module$state.hybrid)[1]+table(gsc.module$state.hybrid)[2])
pheatmap::pheatmap(n,show_colnames = F,annotation_col = colgroup,
                   annotation_colors = anno_col_color,
                   color=colorRampPalette(c("#1E90FF","white","#dc143c"))(100),
                   cluster_cols = F,cluster_rows = F,
                   gaps_col = break_point)


