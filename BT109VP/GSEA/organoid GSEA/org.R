setwd("~/lab/BT109VP/org")
rm(list=ls())
library(Seurat)
library(dplyr)
int.seu <- qs::qread("../int.qs")

int.seu <- subset(int.seu,eGFP<=0)

int.seu <- subset(int.seu,subset = orig.ident%in%c("Eng_DMSO","Eng_VP_L","Eng_VP_S","Org"))

int.seu <- subset(int.seu,subset = celltype=="Malignant",invert=T)

int.seu$treat <- ifelse(int.seu$orig.ident=="Org","org","eng")

############################
int.seu <- subset(int.seu,subset = orig.ident%in%c("Eng_DMSO","Org"))

int.seu <- standard10X(int.seu)

DEG <- FindMarkers(int.seu1,ident.1 = "eng",ident.2 = "org",logfc.threshold = 0)

write.csv(DEG,file = "DEG.ALL.csv")

int.seu$cell_orig <- paste0(int.seu$orig.ident,"_",int.seu$celltype)
cellfordeg<-levels(int.seu$celltype)
Idents(int.seu) <- int.seu$cell_orig
DEG_RESULT <- NULL
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(int.seu, ident.1 = paste0("Eng_VP_L_",cellfordeg[i]), 
                         ident.2 = paste0("Eng_DMSO_",cellfordeg[i]), verbose = FALSE,logfc.threshold = 0.5)
  CELLDEG$celltype <- cellfordeg[i]
  CELLDEG$symbol <- rownames(CELLDEG)
  DEG_RESULT <- rbind(DEG_RESULT,CELLDEG)
}

for (i in unique(DEG_RESULT$celltype)){
  DEG <- DEG_RESULT %>% filter(celltype==i)
  ##Creating Gene List for GSEA
  genelist <- DEG %>% arrange(desc(avg_log2FC))
  gsea_input <- genelist$avg_log2FC
  names(gsea_input) <- genelist$symbol
  ##Select your Gene Set
  dir='/home/hzg/lab/Bulk_analysis/MsigDB/symbol'
  gmts <- list.files(dir,pattern = 'gmt')
  gmts
  #Start GSEA Analysis
  gsea_results <- lapply(gmts, function(gmtfile){
    # gmtfile=gmts[2]
    geneset <- read.gmt(file.path(dir,gmtfile)) 
    print(paste0('Now process the ',gmtfile))
    egmt <- GSEA(gsea_input, TERM2GENE=geneset,eps = 0, verbose=FALSE)
    head(egmt)
    return(egmt)
  })
  # gsea_results[[1]] <- setReadable(gsea_results[[1]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  # gsea_results[[2]] <- setReadable(gsea_results[[2]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
  write.csv(gsea_results_df,file = paste0("./GSEA/GSEA_RESULT_",i,".csv") )
  save(gsea_results,gsea_results_df,file = paste0("./GSEA/GSEA_RESULT_",i,".Rdata"))
}
