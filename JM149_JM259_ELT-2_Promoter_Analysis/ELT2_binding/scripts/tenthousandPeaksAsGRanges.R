library(GenomicRanges)
peaks_df = read.table('allStagesUNION.IDR_0.05.sorted.bed_s.df', header=T)
peaks = makeGRangesFromDataFrame(peaks_df, starts.in.df.are.0based=T,keep.extra.columns=T)
Seqinfo(
c('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrX'),
c(15072434, 15279421, 13783801, 17493829, 20924180, 17718942),
c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
c("ce11", "ce11", "ce11", "ce11", "ce11", "ce11")) -> ce11_info
seqinfo(peaks) <- ce11_info
