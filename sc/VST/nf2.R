
tbl <- xtabs(Freq ~ NF2_mutation + celltype, vst3@meta.data)
proportions(tbl, "Gender")
seurat_colors <- sample(MetBrewer::met.brewer("Klimt",9),9)
cell.prop <- as.data.frame(prop.table(table(vst3$celltype)))
colnames(cell.prop)<- c("sample","proportion")
label_1 <- paste0(substr(cell.prop$proportion*100,1,4),"%")
ggplot(cell.prop,aes(sample,proportion,fill=sample))+
  scale_fill_manual(values = seurat_colors)+
  geom_bar(stat="identity")+
  geom_text(aes(label=label_1))+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))
# int.seu$celltype <- gsub("Myeloid cells","Myeloid",int.seu$celltype)
# int.seu$celltype <- gsub("Schwannoma cells","Schwanncells",int.seu$celltype)
# int.seu$celltype <- gsub("myeSC","Schwanncells",int.seu$celltype)
# int.seu$celltype <- gsub("nmSC","Schwanncells",int.seu$celltype)
DT2 <-vst3@meta.data %>% group_by(NF2_status) %>% tally() %>% as.data.frame()
DT <-vst3@meta.data %>% group_by(NF2_status,celltype) %>% tally()%>% as.data.frame()
DT[1:9,3] <- DT[1:9,3]/DT2[1,2]
DT[10:18,3] <- DT[10:18,3]/DT2[2,2]
DT[19:27,3] <- DT[19:27,3]/DT2[3,2]
DT[28:36,3] <- DT[28:36,3]/DT2[4,2]

cell.prop <- DT
label_1 <- paste0(substr(cell.prop$n*100,1,4),"%")
ggplot(cell.prop,aes(celltype,n,fill=celltype))+
  scale_fill_manual(values = seurat_colors)+
  geom_bar(stat="identity")+
  geom_text(aes(celltype,n+0.02,label=label_1))+
  theme_bw()+theme(axis.ticks.length=unit(0.1,'cm'),
                   title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  guides(fill=guide_legend(title=NULL))+
  facet_wrap(vars(NF2_status))
