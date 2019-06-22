if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (!requireNamespace("DiffBind", quietly = TRUE)) {
    BiocManager::install("DiffBind")
}
library(DiffBind)

if (!requireNamespace("pheatmap", quietly = TRUE)) {
    BiocManager::install("pheatmap")
}
library(pheatmap)

threshold=0.05
# normalized by differential analysis
bNormed=T

if (! "VD" %in% ls()) {
    VD=dba.load('VD',pre='')
}

# contrast 1: LE versus L1
counts_1 = dba.report(VD, contrast = 1, th=threshold, bCalled=T, bNormalized=FALSE)
normed_1 = dba.report(VD, contrast = 1, th=threshold, bCalled=T, bNormalized=bNormed)
# assign column names that are meaningful in all comparisons
counts_1$Fold_LE_L1 = counts_1$Fold
counts_1$Fold = NULL
counts_1$FDR_LE_L1 = counts_1$FDR
counts_1$FDR = NULL
counts_1$Called_LE = counts_1$Called1
counts_1$Called1 = NULL
counts_1$Called_L1 = counts_1$Called2
counts_1$Called2 = NULL
embryo_factor_1 = counts_1$Conc_Embryo - normed_1$Conc_Embryo
larval_1_factor_1 = counts_1$Conc_Larval_1 - normed_1$Conc_Larval_1

# contrast 2: LE versus L3
counts_2 = dba.report(VD, contrast = 2, th=threshold, bCalled=T, bNormalized=FALSE)
normed_2 = dba.report(VD, contrast = 2, th=threshold, bCalled=T, bNormalized=bNormed)
counts_2$Fold_LE_L3 = counts_2$Fold
counts_2$Fold = NULL
counts_2$FDR_LE_L3 = counts_2$FDR
counts_2$FDR = NULL
counts_2$Called_LE = counts_2$Called1
counts_2$Called1 = NULL
counts_2$Called_L3 = counts_2$Called2
counts_2$Called2 = NULL
embryo_factor_2 = counts_2$Conc_Embryo - normed_2$Conc_Embryo
larval_3_factor_2 = counts_2$Conc_Larval_3 - normed_2$Conc_Larval_3

# contrast 3: L1 versus L3
counts_3 = dba.report(VD, contrast = 3, th=threshold, bCalled=T, bNormalized=FALSE)
normed_3 = dba.report(VD, contrast = 3, th=threshold, bCalled=T, bNormalized=bNormed)
counts_3$Fold_L1_L3 = counts_3$Fold
counts_3$Fold = NULL
counts_3$FDR_L1_L3 = counts_3$FDR
counts_3$FDR = NULL
counts_3$Called_L1 = counts_3$Called1
counts_3$Called1 = NULL
counts_3$Called_L3 = counts_3$Called2
counts_3$Called2 = NULL

larval_1_factor_3 = counts_3$Conc_Larval_3 - normed_3$Conc_Larval_3
larval_3_factor_3 = counts_3$Conc_Larval_3 - normed_3$Conc_Larval_3

# take the means of the pairwise normalization factors
embryo_factor = mean(c(embryo_factor_1, embryo_factor_2))
larval_1_factor = mean(c(larval_1_factor_1, larval_1_factor_3))
larval_3_factor = mean(c(larval_3_factor_2, larval_3_factor_3))

# drop unused columns
# p-value
counts_1$`p-value` = NULL
counts_2$`p-value` = NULL
counts_3$`p-value` = NULL
normed_1$`p-value` = NULL
normed_2$`p-value` = NULL
normed_3$`p-value` = NULL
# Conc
counts_1$`Conc` = NULL
counts_2$`Conc` = NULL
counts_3$`Conc` = NULL
normed_1$`Conc` = NULL
normed_2$`Conc` = NULL
normed_3$`Conc` = NULL

# merge
(counts_all = merge(counts_1, counts_2, counts_3, all=TRUE))

peaks = counts_all[, c('FDR_LE_L1', 'FDR_LE_L3','FDR_L1_L3','Fold_LE_L1', 'Fold_LE_L3','Fold_L1_L3', 'Conc_Embryo','Conc_Larval_1','Conc_Larval_3', 'Called_LE','Called_L1','Called_L3')]

# order correctly by chrmosome
seqlevels(peaks) <- sort(seqlevels(peaks))
peaks = sort(peaks)

# apply the normalization factors by subtraction
peaks$Conc_Embryo = peaks$Conc_Embryo - embryo_factor
peaks$Conc_Larval_1 = peaks$Conc_Larval_1 - larval_1_factor
peaks$Conc_Larval_3 = peaks$Conc_Larval_3 - larval_3_factor
# update the Fold changes by subtracting the 'Conc' columns
peaks$Fold_LE_L3 = peaks$Conc_Embryo - peaks$Conc_Larval_3
peaks$Fold_L1_L3 = peaks$Conc_Larval_1 - peaks$Conc_Larval_3
peaks$Fold_LE_L1 = peaks$Conc_Embryo - peaks$Conc_Larval_1
# Consider FDRs of NA to be 1
peaks$FDR_LE_L1[ is.na(peaks$FDR_LE_L1)] <- 1
peaks$FDR_LE_L3[ is.na(peaks$FDR_LE_L3)] <- 1
peaks$FDR_L1_L3[ is.na(peaks$FDR_L1_L3)] <- 1

# a class of Embryo Specific/High
embryo_peaks = peaks[peaks$FDR_LE_L3 < 0.01 & peaks$FDR_LE_L1 < 0.01 & peaks$Fold_LE_L1 > 0 & peaks$Fold_LE_L3 > 0]


peaks_df = as.data.frame(mcols(peaks[, c('Conc_Embryo','Conc_Larval_1','Conc_Larval_3')]))
peaks_df[is.na(peaks_df)] <- 0
#phk = pheatmap(peaks_df, cluster_rows=T, cluster_cols = F, scale="row", show_rownames=F)

# need to check for existence of the file
valerie_peaks_processed_0.050 = read.table('valerie_peaks_processed_0.050.bedlike', header=T, sep = "\t", comment.char='') 
vpp_05_gr = makeGRangesFromDataFrame(valerie_peaks_processed_0.050, keep.extra.columns = T, seqnames.field = 'X.chrom', start.field = 'chromStart', end.field='chromEnd', starts.in.df.are.0based = T)
hits = findOverlaps(peaks, vpp_05_gr)
peaks_in_gr = peaks[from(hits)]
vpp_05_gr_in_peaks = vpp_05_gr[to(hits)]
peaks_in_gr$other_range = ranges(vpp_05_gr_in_peaks)
mcols(peaks_in_gr) <- cbind(mcols(peaks_in_gr), mcols(vpp_05_gr_in_peaks))

hit_vec = findOverlaps(peaks, vpp_05_gr, select="arbitrary")
peaks_not_overlapping_ix = is.na(hit_vec)
peaks_not_overlapping = peaks[peaks_not_overlapping_ix]

sum_L1_1 = peaks_in_gr$mean_log_L1_1*peaks_in_gr$N_log_L1_1
sum_L1_2 = peaks_in_gr$mean_log_L1_2*peaks_in_gr$N_log_L1_2
sum_LE_1 = peaks_in_gr$mean_log_LE_1*peaks_in_gr$N_log_LE_1
sum_LE_2 = peaks_in_gr$mean_log_LE_2*peaks_in_gr$N_log_LE_2
sum_L3_1 = peaks_in_gr$mean_log_L3_1*peaks_in_gr$N_log_L3_1
sum_L3_2 = peaks_in_gr$mean_log_L3_2*peaks_in_gr$N_log_L3_2

sum_LE = apply(cbind(sum_LE_1, sum_LE_2), 1, mean)
sum_L1 = apply(cbind(sum_L1_1, sum_L1_2), 1, mean)
sum_L3 = apply(cbind(sum_L3_1, sum_L3_2), 1, mean)
mean_LE = apply(cbind(peaks_in_gr$mean_log_LE_1,peaks_in_gr$mean_log_LE_2), 1, mean)
mean_L1 = apply(cbind(peaks_in_gr$mean_log_L1_1,peaks_in_gr$mean_log_L1_2), 1, mean)
mean_L3 = apply(cbind(peaks_in_gr$mean_log_L3_1,peaks_in_gr$mean_log_L3_2), 1, mean)
max_LE = apply(cbind(peaks_in_gr$max_log_LE_1,peaks_in_gr$max_log_LE_2), 1, mean)
max_L1 = apply(cbind(peaks_in_gr$max_log_L1_1,peaks_in_gr$max_log_L1_2), 1, mean)
max_L3 = apply(cbind(peaks_in_gr$max_log_L3_1,peaks_in_gr$max_log_L3_2), 1, mean)

# correlation with the Conc_... vectors is higher with sum than max
cor(data.frame(sum_LE, mean_LE, max_LE, peaks_in_gr$Conc_Embryo),use="complete.obs")
cor(data.frame(sum_L1, mean_L1, max_L1, peaks_in_gr$Conc_Larval_1),use="complete.obs")
cor(data.frame(sum_L3, mean_L3, max_L3, peaks_in_gr$Conc_Larval_3),use="complete.obs")

peaks_in_gr$sum_LE = sum_LE
peaks_in_gr$sum_L1 = sum_L1
peaks_in_gr$sum_L3 = sum_L3
peaks_in_gr$max_LE = max_LE
peaks_in_gr$max_L1 = max_L1
peaks_in_gr$max_L3 = max_L3

THRESHOLD = .2
action = peaks_in_gr[peaks_in_gr$peaks_quantile >= THRESHOLD, c('max_LE','max_L1','max_L3', 'kclust_mapping', 'Conc_Embryo', 'Conc_Larval_1', 'Conc_Larval_3')]
action$kclust_mapping[is.na(action$kclust_mapping) ] <- 0
Embryo_min = min(0, min(action$Conc_Embryo, na.rm = T))
Larval_1_min = min(0, min(action$Conc_Larval_1, na.rm = T))
Larval_3_min = min(0, min(action$Conc_Larval_3, na.rm = T))
action$Conc_Embryo[is.na(action$Conc_Embryo)] <- Embryo_min
action$Conc_Larval_1[is.na(action$Conc_Larval_1)] <- Larval_1_min
action$Conc_Larval_3[is.na(action$Conc_Larval_3)] <- Larval_3_min

# add secondary clusters
action$kclust_mapping_2 = 0
#k = 0
#action$kclust_mapping_2[ action$kclust_mapping == k] = kmeans(mcols(action[action$kclust_mapping==k,c('Conc_Embryo', 'Conc_Larval_1', 'Conc_Larval_3')]),4)$cluster
k = 1
action$kclust_mapping_2[ action$kclust_mapping == k] = kmeans(mcols(action[action$kclust_mapping==k,c('Conc_Embryo', 'Conc_Larval_1', 'Conc_Larval_3')]),4)$cluster
k = 2
action$kclust_mapping_2[ action$kclust_mapping == k] = kmeans(mcols(action[action$kclust_mapping==k,c('Conc_Embryo', 'Conc_Larval_1', 'Conc_Larval_3')]),4)$cluster
k = 3
action$kclust_mapping_2[ action$kclust_mapping == k] = kmeans(mcols(action[action$kclust_mapping==k,c('Conc_Embryo', 'Conc_Larval_1', 'Conc_Larval_3')]),4)$cluster
k = 4
action$kclust_mapping_2[ action$kclust_mapping == k] = kmeans(mcols(action[action$kclust_mapping==k,c('Conc_Embryo', 'Conc_Larval_1', 'Conc_Larval_3')]),4)$cluster

# plot
pheatmap(as.data.frame(mcols(action[order(action$kclust_mapping, action$kclust_mapping_2)])), cluster_cols = F,  cluster_rows=F, show_rownames=F)

if (FALSE) {
    stopifnot(file.exists('spp_peaks'))
    source('scripts/readIDR.R');
    LE_IDR
    L1_IDR
    L3_IDR

    # map IDR status onto peaks object
    peaks$IDR_LE = 0
    peaks$IDR_L1 = 0
    peaks$IDR_L3 = 0
    peaks$IDR_LE[ from(findOverlaps(peaks, LE_IDR)) ] = 1
    peaks$IDR_L1[ from(findOverlaps(peaks, L1_IDR)) ] = 1
    peaks$IDR_L3[ from(findOverlaps(peaks, L3_IDR)) ] = 1
    peaks$IDR = apply(cbind(peaks$IDR_LE,peaks$IDR_L1,peaks$IDR_L3), 1,sum)

    prev_df = read.table("allStagesUNION.IDR_0.05.sorted.bed_s.df", header = T, sep="\t")
    prev_gr = makeGRangesFromDataFrame(prev_df, keep.extra.columns = T, starts.in.df.are.0based = T)
}
