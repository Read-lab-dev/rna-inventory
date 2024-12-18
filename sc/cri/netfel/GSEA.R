#####GSEA#####
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(clusterProfiler)
library(GSEABase)
gsc.seu <- readRDS("../gsc.seu.Rds")
Idents(gsc.seu) <- gsc.seu$state

DEG_ALL <- FindAllMarkers(gsc.seu, test.use = "MAST",verbose = FALSE,min.pct = 0.1,
                          logfc.threshold = 0,return.thresh = 1)
colnames(DEG_ALL)[1]<- "symbol"
library(stringr)
DEG_ALL$symbol <- str_split(rownames(DEG_ALL),"[.]",simplify =T)[,1]


###############GSEA#########

for (celltype in levels(gsc.seu)){
  DEG <- DEG_ALL %>% filter(cluster==celltype)
  ##Creating Gene List for GSEA
  genelist <- bitr(DEG$symbol,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(DEG, by= c("SYMBOL"="symbol")) %>% 
    arrange(desc(avg_log2FC))
  
  gsea_input <- genelist$avg_log2FC
  
  names(gsea_input) <- genelist$ENTREZID
  
  ##Select your Gene Set
  
  dir='/home/hzg/rna/Bulk_analysis/MsigDB/'
  
  gmts <- list.files(dir,pattern = 'gmt')
  
  gmts
  
  #Start GSEA Analysis
  
  gsea_results <- lapply(gmts, function(gmtfile){
    # gmtfile=gmts[2]
    geneset <- read.gmt(file.path(dir,gmtfile)) 
    print(paste0('Now process the ',gmtfile))
    egmt <- GSEA(gsea_input, TERM2GENE=geneset, verbose=FALSE)
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
  save(gsea_results,gsea_results_df,file = paste0("GSEA_RESULT_",celltype,".Rdata"))
}



#############Fig7C#########
gsea_kegg_list<-gsea_results_list[[1]]
gsea_hallmark_list<-gsea_results_list[[4]]
gsea_go_list<-gsea_results_list[[3]]
gsea_reactome_list<-gsea_results_list[[2]]
write.csv(gsea_go_list,file = "gsea_go_list.csv")
write.csv(gsea_reactome_list,file = "gsea_reactome_list.csv")
write.csv(gsea_hallmark_list,file = "gsea_hallmark_list.csv")
write.csv(gsea_kegg_list,file = "gsea_kegg_list.csv")
gsea_go<-gsea_results[[3]]
gsea_reactome<-gsea_results[[3]]
gsea_kegg<-gsea_results[[1]]
library(ggsci)
library(ggplot2)
library(ggthemes)
library(clusterProfiler)
library(RColorBrewer)
library(MetBrewer)
ys<-rep(pal_lancet(palette = c("lanonc"),alpha = 0.7)(9),3)

ys <-rep(brewer.pal(12, "Set3"),2)
ys <- rep(met.brewer("Cross",n=7,type="discrete"),3)

p1<-ridgeplot(gsea_c2,showCategory = 20,core_enrichment=T,fill = "ID")+theme_few()+
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        legend.position = "none")+scale_fill_manual(values = ys)+
  scale_y_discrete(labels=function(x)substr(x,10,60))+
  labs(title = "REACTOME")
p1
ggsave(p1,filename = "REACTOME.pdf",height = 10,width = 7)

p2<-ridgeplot(gsea_go,showCategory = 20,core_enrichment=T,fill = "ID")+theme_few()+
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        legend.position = "none")+scale_fill_manual(values = ys)+
  scale_y_discrete(position = "right",labels=function(x)substr(x,3,53))+
  labs(title = "GO")
p2
ggsave(p2,filename = "GO.pdf",height = 10,width = 7)

p3<-ridgeplot(gsea_kegg,showCategory = 20,core_enrichment=T,fill = "ID")+theme_few()+
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        legend.position = "none")+scale_fill_manual(values = ys)+
  scale_y_discrete(position = "right",labels=function(x)substr(x,6,50))+
  labs(title = "KEGG")

p3
ggsave(p3,filename = "KEGG.pdf",height = 10,width = 7)

categorys <- c("KEGG_MAPK_SIGNALING_PATHWAY", "KEGG_VEGF_SIGNALING_PATHWAY",
               "HALLMARK_HEDGEHOG_SIGNALING", "HALLMARK_G2M_CHECKPOINT","HALLMARK_OXIDATIVE_PHOSPHORYLATION",
               "REACTOME_PI3K_AKT_SIGNALING_IN_CANCER","GOBP_CELL_CELL_SIGNALING_BY_WNT","REACTOME_SIGNALING_BY_GPCR",
               "REACTOME_TCR_SIGNALING","REACTOME_SIGNALING_BY_NOTCH4","KEGG_CALCIUM_SIGNALING_PATHWAY",
               "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53")
gsea_results_df$ngene<- as.numeric(substr(gsea_results_df$leading_edge,6,7))*gsea_results_df$setSize*0.01

gsea_plot_data <- gsea_results_df[categorys,] %>% arrange(NES)

gsea_plot_data$name <- factor(gsea_plot_data$ID, levels=unique(gsea_plot_data$ID))

ggplot(data = gsea_plot_data,aes(x=name,y=setSize))+
  geom_bar(aes(alpha=0.7,width=0.6),stat = 'identity')+
  geom_bar(aes(y=ngene,fill= NES,alpha=0.8,width=0.6),stat = 'identity')+
  scale_fill_gradientn(colors = met.brewer("Benedictus"))+theme_bw()+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(hjust = 1))+
  coord_flip()+xlab(NULL)+ylab(NULL)
ggsave("gsea_barplot.pdf")


############Fig7D#########
library(enrichplot)
library(ggsci)

gsea_c2 <- gsea_results[[1]]

gsea_go <- gsea_results[[2]]

ys <- rep(pal_npg("nrc", alpha = 1)(10),6)
# gseaplot2(gsea_go,geneSetID = c("GOBP_GLIAL_CELL_DEVELOPMENT",
#                                 "GOBP_GLIOGENESIS",
#                                 "GOBP_GLIAL_CELL_DIFFERENTIATION"),pvalue_table = F,
#           color = ys,base_size = 9, rel_heights = c(2, 0.5, 1))

gseaplot2(gsea_c2,geneSetID = c("WP_FERROPTOSIS",
                                "HALLMARK_P53_PATHWAY"),pvalue_table = F,
          color = ys,base_size = 9, rel_heights = c(2, 0.5, 1))

gseaplot2(gsea_c2,geneSetID = c("REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES"),pvalue_table = F,
          color = ys,base_size = 9, rel_heights = c(2, 0.5, 1))
gseaplot2(gsea_c2,geneSetID = 'HALLMARK_HYPOXIA',
          title = 'HALLMARK_HYPOXIA', color = "red",
          base_size = 6,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#结果由三张图拼接而成，这里选择展示全部
          pvalue_table = FALSE,#p值表格
          ES_geom = "line")#line or dot


ggsave(filename = paste0("GSEA_",celltype,".pdf"),height = 3,width = 3)

HALLMARK_G2M_CHECKPOINT

gseaplot2(gsea_c2,geneSetID = c("HALLMARK_G2M_CHECKPOINT"),pvalue_table = F,
          color = ys,base_size = 9, rel_heights = c(2, 0.5, 1))
gseaplot2(gsea_go,geneSetID = c("GOBP_GLIAL_CELL_DIFFERENTIATION",
                                "GOBP_GLIOGENESIS"),pvalue_table = F,
          color = ys,base_size = 9, rel_heights = c(2, 0.5, 1))


#############GSEA_GSC###########
gseaplot2(gsea_c2,geneSetID = c("HALLMARK_HEDGEHOG_SIGNALING"),
          title = "HALLMARK_HEDGEHOG_SIGNALING",pvalue_table = F,
          color = "navyblue",base_size = 6, rel_heights = c(2, 0.5, 1))

ggsave("GSEA_GSC_1.pdf",height = 3,width = 3)

gseaplot2(gsea_c2,geneSetID = c("KEGG_MAPK_SIGNALING_PATHWAY"),
          title = "KEGG_MAPK_SIGNALING_PATHWAY",pvalue_table = F,
          color = "navyblue",base_size = 6, rel_heights = c(2, 0.5, 1))

ggsave("GSEA_GSC_2.pdf",height = 3,width = 3)

gseaplot2(gsea_c2,geneSetID = c("HALLMARK_G2M_CHECKPOINT"),
          title = "HALLMARK_G2M_CHECKPOINT",pvalue_table = F,
          color = "firebrick",base_size = 6, rel_heights = c(2, 0.5, 1))
ggsave("GSEA_GSC_3.pdf",height = 3,width = 3)
