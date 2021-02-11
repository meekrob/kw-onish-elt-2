library(GenomicRanges, warn.conflicts=F)
library(dplyr, warn.conflicts=F)
match_ranges = function(dfG, LE_IDR, L1_IDR, L3_IDR) {
  LE_verify = as.data.frame(findOverlaps(LE_IDR, dfG, type="any"))
  L1_verify = as.data.frame(findOverlaps(L1_IDR, dfG, type="any"))
  L3_verify = as.data.frame(findOverlaps(L3_IDR, dfG, type="any"))
  colnames(LE_verify) <- c('index_of_LE_IDR','index_of_dfG')
  colnames(L1_verify) <- c('index_of_L1_IDR','index_of_dfG') 
  colnames(L3_verify) <- c('index_of_L3_IDR','index_of_dfG') 
  # create a full data frame where the IDR hits i's are matched to the merged i, 
  # introduces NAs where there is none in a given set.
  # duplicate rows for some cases where there is a many:one relationship
  IDR_match_i =full_join(L1_verify,
                         full_join(LE_verify,
                                   L3_verify, by="index_of_dfG") ,by="index_of_dfG")
  #IDR_match_i = IDR_match_i[order(IDR_match_i$index_of_dfG),]
  IDR_match_i = IDR_match_i[,c('index_of_dfG','index_of_LE_IDR','index_of_L1_IDR','index_of_L3_IDR')]
  return(IDR_match_i)
}

score_peak_summit_agreement = function(dfG, LE_IDR, L1_IDR, L3_IDR) {
  IDR_match_i = match_ranges(dfG, LE_IDR, L1_IDR, L3_IDR)
  LE_IDR_summits = data.frame(
    index_in_dfG = to(LE_verify),
    LE_peak_1_summit = LE_IDR$peak_1_start + LE_IDR$peak_1_summit,
    LE_peak_2_summit = LE_IDR$peak_2_start + LE_IDR$peak_2_summit
  )

  L1_IDR_summits = data.frame(
    index_in_dfG = to(L1_verify),
    L1_peak_1_summit = L1_IDR$peak_1_start + L1_IDR$peak_1_summit,
    L1_peak_2_summit = L1_IDR$peak_2_start + L1_IDR$peak_2_summit
  )
  L3_IDR_summits = data.frame(
    index_in_dfG = to(L3_verify),
    L3_peak_1_summit = L3_IDR$peak_1_start + L3_IDR$peak_1_summit,
    L3_peak_2_summit = L3_IDR$peak_2_start + L3_IDR$peak_2_summit
  )
  
  IDR_summits = full_join(full_join(LE_IDR_summits,L1_IDR_summits,by='index_in_dfG'), L3_IDR_summits, by='index_in_dfG')
  #IDR_summits = IDR_summits[order(IDR_summits$index_in_dfG),]
  # average number of bp any summit is from the overall ones in a IDR peak
  summit_mean = function(x) {sum(abs(x - mean(x, na.rm=T)),na.rm=T)/sum(!is.na(x))}
  # summit_score = function(x) { (6/sum(!is.na(x)))*summit_mean(x)} # scale inversely by missing data
  stds=apply(IDR_summits[,2:7], 1, summit_mean)
  reduced_mapping = match(1:length(dfG), IDR_summits$index_in_dfG) # if there are many-to-one mappings, the first encountered is used
  summit_agreement = stds[ reduced_mapping ] 
  IDR_summits_reduced = IDR_summits[ reduced_mapping,]
  IDR_summits_reduced$summit_agreement = summit_agreement
  return(IDR_summits_reduced)
}

jaccard = function(LE_IDR,L1_IDR,L3_IDR,LE_IDR_i,L1_IDR_i,L3_IDR_i){
  # go through the merged dataset, one row at a time, and get the ranges that comprise the interval for that row
  jaccard = c()
  for (i_dfG in 1:length(dfG)) {
    LE_IDR_match_i = IDR_match_i[IDR_match_i$index_of_dfG==i_dfG,2]
    L1_IDR_match_i = IDR_match_i[IDR_match_i$index_of_dfG==i_dfG,3]
    L3_IDR_match_i = IDR_match_i[IDR_match_i$index_of_dfG==i_dfG,4]
    LE_IDR_match_i = unique(LE_IDR_match_i[!is.na(LE_IDR_match_i)])
    L1_IDR_match_i = unique(L1_IDR_match_i[!is.na(L1_IDR_match_i)])
    L3_IDR_match_i = unique(L3_IDR_match_i[!is.na(L3_IDR_match_i)])
    ranges=c(LE_IDR[LE_IDR_match_i], L1_IDR[L1_IDR_match_i],L3_IDR[L3_IDR_match_i])
    tryCatch({
      isect_width = GenomicRanges::width(r_intersect(ranges))
      union_width = GenomicRanges::width(r_union(ranges))
      jaccard = c(jaccard, isect_width/union_width)
    }, error = function(e) { print(i_dfG); jaccard = c(jaccard,0)})
  } 
  return(jaccard)
}

r_intersect = function(ranges) {
  n = length(ranges)
  if (n > 2) {
    return(GenomicRanges::intersect(ranges[1], r_intersect(ranges[2:n])));
  }
  else {
    return(GenomicRanges::intersect(ranges[1],ranges[2]))
  }
}

r_union = function(ranges) {
  n = length(ranges)
  if (n > 2) {
    return(GenomicRanges::union(ranges[1], r_union(ranges[2:n])));
  }
  else {
    return(GenomicRanges::union(ranges[1],ranges[2]))
  }
}