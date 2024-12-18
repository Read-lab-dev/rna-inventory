
library(dplyr)
library(clusterProfiler)
 
for (i in unique(dmso.deg$cluster)){
  dat <- dmso.deg %>% filter(cluster==i)
  dat <- dat %>% arrange(desc(avg_log2FC))
  gsea_input <- dat$avg_log2FC
  names(gsea_input) <- dat$gene
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
  write.csv(gsea_results_df,file = paste0("./GSEA_RESULT1_",i,".csv") )
  save(gsea_results,gsea_results_df,file = paste0("./GSEA_RESULT1_",i,".Rdata"))
}
