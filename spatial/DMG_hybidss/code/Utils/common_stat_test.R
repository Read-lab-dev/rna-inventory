# fisher_test <- function(a, b, bg){
#   a = unlist(a); b = unlist(b)
#   x = length(intersect(a,b))
#   m = length(setdiff(a,b))
#   n = length(setdiff(b,a))
#   k = length(setdiff(bg, union(a,b)))
#   mat = matrix(c(x,m,n,k),2,2)
#   print(mat)
#   return(fisher.test(mat)$p.value)
# }

## Fisher exact test to test over-representaiton of gene sets
## @param a: genes of interest, e.g. metaprogram
## @param b: gene program to compare, e.g. GO term 
## @param bg: background genes (e.g. all genes)
fisher_test <- function(a, b, bg){
  a = unlist(a); b = unlist(b)
  x = length(intersect(a, b))
  m = length(a) - x
  n = length(intersect(bg, b)) - x
  k = length(bg) - n - length(a)
  mat = matrix(c(x,m,n,k),2,2)
  ##print(mat)
  pval = fisher.test(mat, alternative = "greater")$p.value
  res = c(mat[1, 1], colSums(mat)[1], mat[1, 1]/colSums(mat)[1],
          mat[1, 2], colSums(mat)[2], mat[1, 2]/colSums(mat)[2],
          pval)
  names(res) = c("deg_hit", "deg_total", "deg_ratio",
                 "bg_hit", "bg_total", "bg_ratio", "pvalue")
  return(res)
}

jaccard_index <- function(a, b){
    a = unlist(a); b=unlist(b)
    return(length(intersect(a,b))/length(union(a,b)))
}