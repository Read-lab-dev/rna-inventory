computeEucDist <- function(sample_name){
  mat = NES_all_samples[[sample_name]]
  mat_dist = dist(t(mat))
}

calinsky <- function(hhc, dist = NULL, gMax = round(1 + 3.3 * log(length(hhc$order), 10))) {
  msg <- ""
  if (is.null(dist)) {
    require(clue)
    dist <- sqrt(as.cl_ultrametric(hhc))
    # message(msg <- "The distance matrix not is provided, using the cophenetic matrix")
  } else if (attr(dist, "method") != "euclidean") {
    require(clue)
    dist <- sqrt(as.cl_ultrametric(hhc))
    # message(msg <- "The distance matrix is not euclidean, using the cophenetic matrix")
  }
  dist <- as.matrix(dist)^2
  A <- -dist/2
  A_bar <- apply(A, 1, mean)
  totalSum <- sum(diag(A) - 2 * A_bar + mean(A))
  n <- length(hhc$order)
  ans <- rep(0, gMax)
  for (g in 2:gMax) {
    cclust <- cutree(hhc, k = g)
    withinSum <- 0
    for (k in 1:g) {
      if (sum(cclust == k) == 1) 
        next
      A <- as.matrix(-dist/2)[cclust == k, cclust == k]
      A_bar <- apply(A, 1, mean)
      withinSum <- withinSum + sum(diag(A) - 2 * A_bar + 
                                     mean(A))
    }
    betweenSum <- totalSum - withinSum
    betweenSum <- betweenSum/(g - 1)
    withinSum <- withinSum/(n - g)
    ans[g] <- betweenSum/withinSum
  }
  class(ans) <- "calinsky"
  attr(ans, "message") <- msg
  return(ans)
}

runConsensusClustering <- function(sample_name, maxK=6, reps=10000, pItem=0.7){
  message(paste0("Processing sample: ", sample_name, 
                 "; Starting time: ", Sys.time(), "\n"))
  print(sample_name)
  tmp = EucDist_all_samples[[sample_name]]
  title = paste0(seurat_fig_folder, sample_name, "_rep", reps)
  cc_res = ConsensusClusterPlus(tmp,
                                maxK=maxK,
                                reps=reps,
                                pItem=pItem,
                                pFeature=1,
                                innerLinkage="ward.D2", 
                                finalLinkage="ward.D2",
                                title=title,
                                clusterAlg="hc",
                                distance="euclidean",
                                seed=123456,plot="png")
  icl = calcICL(cc_res,title=title,plot="png")
  res = list("cc_res" = cc_res, "icl" = icl)
  message("===============================")
  return(res)
}

computeSI <- function(sample_name, clust_index){
  ## Compute dist matrix from consensusClust results
  aConsensus = cc_res_list[[sample_name]][["cc_res"]][[clust_index]]
  consMatrix <- 1 - aConsensus$consensusMatrix
  rownames(consMatrix) <- colnames(consMatrix) <- names(aConsensus$consensusClass) 
  consMatrix <- as.dist(consMatrix)
  res2 = aConsensus$consensusClass
  ## Compute silhouette score
  si2 <- silhouette(x = res2, dist = consMatrix)
  si2 = si2[1:length(res2), ]
  rownames(si2) = names(res2)
  si2 = data.frame(si2)
  return(si2)
}

filterCCres <- function(sample_name, minCell = 10, sil_cutoff = 0.5){
  df = cc_sl_res_list[[sample_name]]
  ## Filter by # of cells in each cluster; remove those in cluster with # of cell < 10
  nCells = table(df$cluster)
  clust_to_keep = names(nCells)[nCells >= minCell]
  df = df[((df$cluster %in% clust_to_keep) & (df$sil_width>=sil_cutoff)), ]
  return(df)
}

runGeneMWW <- function(gene, clust, df_cm, df_clust){
  a = df_cm[gene, df_clust$cluster == clust]
  b = df_cm[gene, df_clust$cluster != clust]
  res = wilcox.test(a, b)
  res2 = res$statistic/(length(a)*length(b))
  res2 = log2(res2/(1-res2))
  return(res2)
}

runSampleGeneMWW <- function(sample_name){
  message(paste0("Processing sample: ", sample, 
                 "; Starting time: ", Sys.time(), "\n"))
  counts = cm_list[[sample_name]]
  df = cc_sl_res_filtered_list[[sample_name]]
  counts = counts[, rownames(df)]
  all_clusters = sort(unique(df$cluster))
  res = NULL
  res_rownames = NULL
  for (gene in rownames(counts)){
    res1 = NULL
    for (clust in all_clusters){
      tmp = runGeneMWW(gene=gene, clust=clust, df_cm=counts, df_clust=df)
      res1 = c(res1, tmp)
    }
    res = rbind(res, res1)
    res_rownames = c(res_rownames, gene)
  }
  rownames(res) = res_rownames
  colnames(res) = all_clusters
  message("=============Done=============")
  return(res)
}

runBinarization <- function(p, fdr, thres=0.58, fdr_thres=0.01){
  res = ifelse((p>=thres) & (fdr<=fdr_thres), 1, 0)
  return(res)
}

computeJD <- function(a, b){
  ## If both clusters have NO active pathways, their distance = 1
  if (length(a) == 0 & length(b) == 0){
    res = 1
  } else{
    res = 1- length(intersect(a, b))/length(union(a, b))
  }
  return(res)
}

findCommonPathways <- function(NES_list){
  pathways = names(NES_list[[1]])
  for (i in NES_list){
    pathways = intersect(pathways, names(i))
  }
  return(pathways)
}

findCommonPathwaysClust <- function(x, cc_res = cc_res_clusters){
  cur_clusters = names(cc_res)[cc_res == x]
  cur_NES_list = NES_clust_bi_list2[cur_clusters]
  common_pathways = findCommonPathways(cur_NES_list)
  cur_NES_list_updated = sapply(cur_NES_list, function(x) x[common_pathways])
  return(cur_NES_list_updated)
}