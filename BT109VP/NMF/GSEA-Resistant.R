library(dplyr)
library(clusterProfiler)
library(Seurat)
dmso.seu <- qs::qread("../gsc.seu.qs")
dmso.seu <- subset(dmso.seu,celltype=="tCYC")
Idents(dmso.seu) <- ifelse(dmso.seu$orig.ident=="Eng_VP_L","resistant","control")
dmso.deg <- FindMarkers(dmso.seu,ident.1 = "resistant",ident.2 = "control",logfc.threshold=0)

  dat <- dmso.deg %>% filter(p_val_adj<=0.05)
  dat <- dat %>% arrange(desc(avg_log2FC))
  gsea_input <- dat$avg_log2FC
  names(gsea_input) <- rownames(dat)
  ##Select your Gene Set
  dir='/home/hzg/rna/Bulk_Analysis/MsigDB/'
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
  # gsea_results[[1]] <- setReadable(gsea_results[[1]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  # gsea_results[[2]] <- setReadable(gsea_results[[2]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
  write.csv(gsea_results_df,file = "./GSEA_RESULT_cyc-Resistant.csv")
  save(gsea_results,gsea_results_df,file ="./GSEA_RESULT_cyc-Resistant.Rdata")
