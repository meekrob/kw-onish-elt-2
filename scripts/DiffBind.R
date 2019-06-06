if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (!requireNamespace("DiffBind", quietly = TRUE)) {
    BiocManager::install("DiffBind")
}
library(DiffBind)

threshold=0.05
# normalized by differential analysis
bNormed=T

if (! "VD" %in% ls()) {
    dba.load("VD.RData")
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

business = counts_all[, c('FDR_LE_L1', 'FDR_LE_L3','FDR_L1_L3','Fold_LE_L1', 'Fold_LE_L3','Fold_L1_L3', 'Conc_Embryo','Conc_Larval_1','Conc_Larval_3', 'Called_LE','Called_L1','Called_L3')]
# apply the normalization factors by subtraction
business$Conc_Embryo = business$Conc_Embryo - embryo_factor
business$Conc_Larval_1 = business$Conc_Larval_1 - larval_1_factor
business$Conc_Larval_3 = business$Conc_Larval_3 - larval_3_factor
# update the Fold changes by subtracting the 'Conc' columns
business$Fold_LE_L3 = business$Conc_Embryo - business$Conc_Larval_3
business$Fold_L1_L3 = business$Conc_Larval_1 - business$Conc_Larval_3
business$Fold_LE_L1 = business$Conc_Embryo - business$Conc_Larval_1
# Consider FDRs of NA to be 1
business$FDR_LE_L1[ is.na(business$FDR_LE_L1)] <- 1
business$FDR_LE_L3[ is.na(business$FDR_LE_L3)] <- 1
business$FDR_L1_L3[ is.na(business$FDR_L1_L3)] <- 1

# a class of Embryo Specific/High
embryo_business = business[business$FDR_LE_L3 < 0.01 & business$FDR_LE_L1 < 0.01 & business$Fold_LE_L1 > 0 & business$Fold_LE_L3 > 0]

business_df = as.data.frame(mcols(business[, c('Conc_Embryo','Conc_Larval_1','Conc_Larval_3')]))
business_df[is.na(business_df)] <- 0
log_business_df = log(business_df + 1)
phk = pheatmap(log_business_df, cluster_rows=T, cluster_cols = F, scale="row")
