#Anno
int.seu <-qs::qread("tumor_adj.qs")
Idents(int.seu) <- int.seu$celltype
all.markers <- FindAllMarkers(int.seu, only.pos = TRUE, min.pct = 0.25,max.cells.per.ident = 200,
                              logfc.threshold = 0) 
DEG_ALL <-all.markers
###############GSEA#########
for (i in unique(DEG_ALL$cluster)){
  DEG <- DEG_ALL %>% filter(cluster==i)
  ##Creating Gene List for GSEA
  genelist <- bitr(DEG$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db") %>% 
    inner_join(DEG, by= c("SYMBOL"="gene")) %>% 
    arrange(desc(avg_log2FC))
  
  gsea_input <- genelist$avg_log2FC
  
  names(gsea_input) <- genelist$ENTREZID
  
  ##Select your Gene Set
  
  dir='/home/hzg/rna/MsigDB/cell/'
  
  gmts <- list.files(dir,pattern = 'gmt')
  
  gmts
  
  #Start GSEA Analysis
  
  gsea_results <- lapply(gmts, function(gmtfile){
    # gmtfile=gmts[2]
    geneset <- read.gmt(file.path(dir,gmtfile)) 
    print(paste0('Now process the ',gmtfile))
    egmt <- GSEA(gsea_input, TERM2GENE=geneset, verbose=FALSE,scoreType = "pos")
    head(egmt)
    return(egmt)
  })
  
  gsea_results[[1]] <- setReadable(gsea_results[[1]], OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
  
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  
  gsea_results_df <- do.call(rbind, gsea_results_list)
  save(gsea_results,gsea_results_df,file = paste0("./anno/GSEA_RESULT_",i,".Rdata"))
}
####

