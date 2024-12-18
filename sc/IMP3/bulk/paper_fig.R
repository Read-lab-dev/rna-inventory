library(dplyr)
library(tibble)
library(MetBrewer)
library(ggrepel)
library(ggplot2)
library(ggthemes)
library(patchwork)
library("RColorBrewer")
library(pheatmap)
library(sva)
library(DESeq2)
resdata <- read.csv("DEseq2_GBM301_shRIOK2.csv",row.names = 1)


resdata$change= as.factor(ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange) > 0.5,
                                 ifelse(resdata$log2FoldChange > 0.5,"up", "down"), "not"))
resdata$change <- factor(resdata$change,levels = c("up","not","down"))

this_tile <- paste0("GBM301_shRIOK2",
                    '\nUp regulated: ',nrow(na.omit(resdata[resdata$change =='up',])) ,
                    '\nDown regulated: ',nrow(na.omit(resdata[resdata$change =='down',])))
label2show <- c("RIOK2","IGF2BP3","MYC","EGFR","OLIG1","OLIG2","SOX4","TOP2A","UBE2C","MKI67","SLC1A3")
label2show <- c("IGF2BP3","MYC","MKI67","TOP2A","UBE2C","PDGFRA","SOX4","SOX9","SOX11")
ggplot(data=na.omit(resdata), aes(x=log2FoldChange, y=-log10(padj), colour=change)) +
  scale_color_manual(values = c('#FF8889','grey80','skyblue'))+
  geom_point(alpha=1, size=0.5) +
  geom_vline(xintercept = c(-0.5,0.5),lty=2,col="azure4",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="azure4",lwd=0.5)+
  geom_label_repel(data = resdata[label2show,],
                   aes(x=log2FoldChange, y=-log10(padj), colour=change,label=gene))+
  theme_set(theme_set(theme_bw(base_size=10)))+
  xlab("log2 Fold Change") + ylab("-log10 adjp-value") + ggtitle(this_tile)+
  theme_clean()+ylim(c(0,100))+
  theme(plot.title = element_text(size=10,hjust = 0.5),
        plot.title.position = "plot",
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.position = c(1, 0.8),
        legend.justification = c(1, 1))
ggsave(filename = "GBM301_shRIOK2_Volcano_For_paper.pdf",width  = 6,height = 4)

load("~/rna/sc/IMP3/bulk/GSEA_GBM301_sh77.Rdata")
library(enrichplot)
gseaplot2(gsea_results[[1]],"YU_MYC_TARGETS_UP",title = "YU_MYC_TARGET_UP",
          color = "cornflowerblue",base_size = 10,rel_heights = c(1.5, 0.5, 0.5))

load("/home/hzg/rna/sc/IMP3/bulk/GSEA_GBM301_sh75.Rdata")
gseaplot2(gsea_results[[1]],geneSetID =c("YU_MYC_TARGETS_UP","VERHAAK_GLIOBLASTOMA_CLASSICAL"),
          title = "GBM301_shIMP3",
          color = c("darkred","cornflowerblue"),base_size = 8,rel_heights = c(1.5, 0.5, 0.5))

gseaplot2(gsea_results[[3]],geneSetID =c("HALLMARK_P53_PATHWAY"),title = "HALLMARK_P53_PATHWAY",
          color = c("cornflowerblue"),base_size = 10,rel_heights = c(1.5, 0.5, 0.5))
