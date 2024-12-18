#################
### Figure 5a ###
#################
FeaturePlot(seurat_obj, 
            features = c("microglia_score", "macrophage_score")
            #cols = c("blue", "lightgrey", "red")
) & NoAxes() &
  scale_color_gradient2(low = "blue", mid = "lightgrey", high="red")

#################
### Figure 5b ###
#################
mg_colors = c("macrophage" = greens[5], "microglia" = greens[8])
DimPlot(object = seurat_obj, 
        group.by = "classification", 
        label = TRUE, 
        cols = mg_colors,
        pt.size = 3, 
        label.size = 8) + NoAxes()

#################
### Figure 5c ###
#################
FeaturePlot(seurat_obj,
            c("P2RY12","TMEM119", 
              "CD163", "TGFBI"),
            cols = c("lightgrey", "red"),
            pt.size = 1,
            ncol = 2) & NoAxes() & NoLegend()

#################
### Figure 5d ###
#################
tmp = table(seurat_obj$sample)
cond = tmp >= 10
df = table(seurat_obj$sample, seurat_obj$classification)
table(names(tmp) == rownames(df))
df = df[cond, ]
df = df/rowSums(df)
df = data.frame(df)
colnames(df) = c("Sample", "Annot", "Proportion")

## Add age
tmp = table(seurat_obj$sample, seurat_obj$peds)
tmp = tmp[cond, ]
tmp2 = apply(tmp, 2, function(x) names(x)[x>0])
df$age = ifelse(df$Sample %in% tmp2$adult, "adult", "pediatric")
df$Annot = gsub("m", "M", df$Annot)

## Add location
tmp = table(seurat_obj$sample, seurat_obj$location)
tmp = tmp[cond, ]
tmp2 = apply(tmp, 2, function(x) names(x)[x>0])
df$location = ifelse(df$Sample %in% tmp2$pontine, "pontine", 
                     ifelse(df$Sample %in% tmp2$thalamic, "thalamic", "others"))

## Plot
ggplot(df, aes(x=age, y=Proportion, fill=age)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 2) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="black", size=1) +
  scale_fill_manual(values=colors_age) + xlab("") +
  facet_wrap(.~Annot) + 
  theme(strip.text = element_text(size=20),
        strip.background = element_blank(),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=20, hjust=1, angle=45),
        axis.text.y = element_text(size=20),
        legend.position = "None")

#################
### Figure 5e ###
#################
p = VlnPlot(seurat_obj, 
            features = c("OSM"), 
            group.by = "classification", 
            split.by = "peds",
            cols = colors_age,
            sort = T, pt.size = 1, ncol = 1)
p & xlab("") & ylab("Log expression level") &
  scale_color_manual(values = colors_age) &
  theme(plot.title = element_text(size=20, face="plain"),
        axis.text = element_text(size=16),
        axis.text.x = element_text(size=16, angle=45, hjust=1),
        axis.text.y = element_text(size=16))

#################
### Figure 5f ###
#################
p = VlnPlot(seurat_obj, 
            features = "OSMR", 
            group.by = "age",
            ncol = 1,
            pt.size = 0.01, 
            cols = colors_age) & NoLegend()
p & xlab("") & ylab("Log expression level") &
  theme(plot.title = element_text(size=20, face="plain"),
        axis.text = element_text(size=16),
        axis.text.x = element_text(size=16, angle=45, hjust=1),
        axis.text.y = element_text(size=16))

#################
### Figure 5g ###
#################
p = VlnPlot(seurat_obj, 
            features = c("VIM", "ANXA2", "CD44", "CHI3L1"), 
            group.by = "peds",
            cols = colors_age,
            sort = T, pt.size = 0.1, ncol = 4)

p & labs(x="", y="Log expression level") & 
  theme(axis.text.x = element_text(size=20, angle=45, hjust=1),
        axis.text.y = element_text(size=20))

#################
### Figure 5h ###
#################
seuratobject$celltype <- factor(seuratobject$celltype, levels=c("E14.5 Microglia","P7 Microglia", "P60 Microglia", "P7 Myeloid","P60 Myeloid"), ordered=TRUE)

x <- AverageExpression(seuratobject, features=c("Osm","Hbegf","Cd44", "Anxa1","Anxa2","Chi3l1", "Vim"), group.by="celltype", return.seurat=TRUE)
x$celltype <- colnames(x)
x$celltype <- factor(x$celltype, levels=c("E14.5 Microglia","P7 Microglia", "P60 Microglia", "P7 Myeloid","P60 Myeloid"), ordered=TRUE)

DoHeatmap(x, features=c("Osm","Hbegf","Cd44", "Anxa1","Anxa2","Chi3l1", "Vim"), group.by="celltype", draw.lines=F)+scale_fill_gradientn(colors = c("#2166AC", "white", "#B2182B"))