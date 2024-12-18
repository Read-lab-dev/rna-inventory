## Function to derive pseudobulk (aggregated gene expressions for each cell subpopulation)
## Based on classification by NMF
## @para seurat_obj seurat object that stores cm
## @para nmf_score df that stores scores of each metaprogram; seurat_obj and nmf_score should have identical cells
## @para out_dir output directory to store results
## @para aggr_method method to aggregate results, either mean or median 
## @score_diff discard cells with metaprogram score less than this cutoff when aggregating 
makePseudobulk <- function(seurat_obj, nmf_score, 
                           out_dir = seurat_analysis_folder,
                           aggr_method="mean", score_cutoff=1){
    metagene_program_names = sort(unique(nmf_score$signature_1))
    pseudobulk = NULL
    
    for (metagene in metagene_program_names){
        cm = seurat_obj@raw.data[,nmf_score$signature_1 == metagene & 
                                     nmf_score$score_1 >= score_cutoff]
        if (aggr_method == "mean"){
            pseudobulk = cbind(pseudobulk, rowMeans(cm))
        }
        if (aggr_method == "median"){
            pseudobulk = cbind(pseudobulk, rowMedians(cm))
        }
    }
    colnames(pseudobulk) = metagene_program_names
    rownames(pseudobulk) = rownames(seurat_obj@raw.data)
    saveRDS(pseudobulk, file=paste0(out_dir, "pseudobulk_", aggr_method, ".RDS"))
}