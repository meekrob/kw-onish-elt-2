stopifnot( all(c("row_scaled_x","threshold_ix","dfG") %in% ls()) )
# Develop a higher-confidence reliable set:
# The following selects the "excluded" peaks that are called in all stages,
# and puts the others back into the main dataset. It is cleaner to 
# steady class
steady_set_means = colMeans(row_scaled_x[ !threshold_ix, ])
steady_set_means = cbind(steady_set_means[1:2],steady_set_means[3:4],steady_set_means[5:6])
colnames(steady_set_means) <- c('LE','L1','L3')
rownames(steady_set_means) <- c('rep1','rep2')
steady_set_means = colMeans(steady_set_means)

# "excluded" means "constant during development", and needs an additional constraint of having IDR at all points
# the peakfiles retain the component peaks from replicates, but the comprehensive range is used
called_in_all = dfG$LE_IDR & dfG$L1_IDR & dfG$L3_IDR
threshold_ix_idr = !threshold_ix & called_in_all
two_of_three_idr = (dfG$LE_IDR & dfG$L1_IDR) | (dfG$LE_IDR & dfG$L3_IDR) | (dfG$L1_IDR & dfG$L3_IDR)
called_in_two_or_three_not_changing = dfG[! threshold_ix & two_of_three_idr]
exclude_ix = ! threshold_ix & called_in_all
include_ix = ! exclude_ix

include_x = row_scaled_x[ include_ix,] # include_x now re-includes any of the filtered peaks that weren't called by IDR in                                        # all 3 stages. This results in a very small number of "not-changing" peaks.
##############