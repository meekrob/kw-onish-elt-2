if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    BiocManager::install("GenomicRanges")
}
library(GenomicRanges)
IDRNarrowPeakColumns=c("seqname","start","end","name", "score", "strand", "signalValue","ignored1","ignored2", "peakSummitOffset", "minusLogPValue","minusLogQValue", "start1", "end1", "signalVal1", "peakSummitOffset1", "start2","end2", "signalVal2", "peakSummitOffset2")


# NAS must be connected
stopifnot(file_test("-d", "spp_peaks"))
df = read.table('spp_peaks/IDR/LE_1_LE_2.IDR_0.05.narrowPeak',col.names=IDRNarrowPeakColumns)
LE_IDR = makeGRangesFromDataFrame(df, keep.extra.columns = T, starts.in.df.are.0based = T)
df = read.table('spp_peaks/IDR/L1_1_L1_2.IDR_0.05.narrowPeak',col.names=IDRNarrowPeakColumns)
L1_IDR = makeGRangesFromDataFrame(df, keep.extra.columns = T, starts.in.df.are.0based = T)
df = read.table('spp_peaks/IDR/L3_1_L3_2.IDR_0.05.narrowPeak',col.names=IDRNarrowPeakColumns)
L3_IDR = makeGRangesFromDataFrame(df, keep.extra.columns = T, starts.in.df.are.0based = T)
# drop weird columns
LE_IDR$ignored1 = NULL
LE_IDR$ignored2 = NULL
L1_IDR$ignored1 = NULL
L1_IDR$ignored2 = NULL
L3_IDR$ignored1 = NULL
L3_IDR$ignored2 = NULL
# add stage ID
LE_IDR$stage="LE"
L1_IDR$stage="L1"
L3_IDR$stage="L3"
# merge all
ALL_IDR = sort(c(LE_IDR, L1_IDR, L3_IDR))
