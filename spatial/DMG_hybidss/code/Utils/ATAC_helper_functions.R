library(GenomicRanges)

## ==========================================================
## Converting peaks in string format into a grange object
## ==========================================================

## Converting peaks in string format into a grange object
## @para: peaks (chr-start-end)
## @returns: a grange object with extra metadata (peaks, peak_idx)
peaksToGranges <- function(peaks){
  peaks = data.frame(peaks)
  peaks$chr = sapply(peaks$peaks, function(x) unlist(strsplit(x, split="-"))[1])
  peaks$start = sapply(peaks$peaks, function(x) unlist(strsplit(x, split="-"))[2])
  peaks$end = sapply(peaks$peaks, function(x) unlist(strsplit(x, split="-"))[3])
  peaks$peak_idx = 1:nrow(peaks)
  peaks = peaks[, c("chr", "start", "end", "peaks", "peak_idx")]
  peak_gr <- peaks %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}