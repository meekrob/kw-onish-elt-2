if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt")
}
library(biomaRt)
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    BiocManager::install("GenomicRanges")
}
library(GenomicRanges)



getCodingGenes = function(peaks, within_genes_kb = 5){
  names(peaks) <- peaks$name
  
  
  # Get only the protein-coding genes. Ranges are comprehensive across splice variants.
  # Get the specs from https://parasite.wormbase.org/biomart/martview/ 
  if (file.exists('genes_coding.RData')) {
    load('genes_coding.RData')
  }
  else {
    paramart <-
    useMart("parasite_mart",
            dataset = "wbps_gene",
            host = "https://parasite.wormbase.org",
            port = 443)
    genes_coding = getBM(
      mart = paramart,
      filter = c("species_id_1010",
                 "biotype"),
      value = list(species_id_1010 = "caelegprjna13758",
                   biotype = "protein_coding"),
      attributes = c(
        'wbps_gene_id',
        'chromosome_name',
        'start_position',
        'end_position',
        'strand',
        'wormbase_gseq',
        "wormbase_locus"
      )
    )
    save(genes_coding, file="genes_coding.RData")
  }
  # exclude mitochondria
  genes_coding_noMT = genes_coding[genes_coding$chromosome_name != 'MtDNA',]
  
  bedformat = genes_coding_noMT[,c('chromosome_name','start_position','end_position','strand','wbps_gene_id')]
  bedformat$strand = with(bedformat, ifelse(strand==1, '+','-'))
  bedformat = bedformat[with(bedformat, order(chromosome_name, start_position)), ] # WS271 as of this writing
  
  # add 'chr'
  genes_coding_noMT$chromosome_name = paste('chr', genes_coding_noMT$chromosome_name, sep='')
  # strands must be '+/-', but paramart returns 1/-1
  strand_plus_minus = genes_coding_noMT$strand
  strand_plus_minus[ strand_plus_minus == '1'] <- '+'
  strand_plus_minus[ strand_plus_minus == '-1'] <- '-'
  genes_coding_noMT$strand = strand_plus_minus
  all_CDS_genes = makeGRangesFromDataFrame(genes_coding_noMT, 
                                           keep.extra.columns = T, 
                                           seqnames.field='chromosome_name', 
                                           start.field="start_position", 
                                           end.field="end_position")
  
  names(all_CDS_genes) <- all_CDS_genes$wbps_gene_id
  all_CDS_genes$name = ifelse(all_CDS_genes$wormbase_locus == '', all_CDS_genes$wormbase_gseq, all_CDS_genes$wormbase_locus)
  
  # do 'overlapping' as an unambiguous criterion first
  overlapping_ap = annotatePeakInBatch(
    peaks, AnnotationData = all_CDS_genes, bindingRegions = c(-within_genes_kb, within_genes_kb),
    featureType = "Exon",
    PeakLocForDistance = "middle",
    FeatureLocForDistance = "middle",
    output = "overlapping")
  peaks[overlapping_ap$peak[is.na(overlapping_ap$feature) ]] -> no_overlap_peaks
 
  # sort by shortest distance within each peak
  overlapping_ap = overlapping_ap[ order(overlapping_ap$peak,overlapping_ap$shortestDistance)]
  
  unique_indexes = function(vec) {
    return(match(unique(vec), vec))
  }
  
  upstream_of_peak_togene_end_ap = annotatePeakInBatch(
    no_overlap_peaks, AnnotationData = all_CDS_genes, bindingRegions = c(-within_genes_kb, within_genes_kb),
    featureType = "Exon",
    PeakLocForDistance = "start",
    FeatureLocForDistance = "end",
    output = "nearestLocation") 
  
  upstream_of_peak_togene_start_ap = annotatePeakInBatch(
    no_overlap_peaks, AnnotationData = all_CDS_genes, bindingRegions = c(-within_genes_kb, within_genes_kb),
    featureType = "Exon",
    PeakLocForDistance = "start",
    FeatureLocForDistance = "start",
    output = "nearestLocation")  
  
  downstream_of_peak_togene_end_ap = annotatePeakInBatch(
    no_overlap_peaks, AnnotationData = all_CDS_genes, bindingRegions = c(-within_genes_kb, within_genes_kb),
    featureType = "Exon",
    PeakLocForDistance = "end",
    FeatureLocForDistance = "end",
    output = "nearestLocation") 
  
  downstream_of_peak_togene_start_ap = annotatePeakInBatch(
    no_overlap_peaks, AnnotationData = all_CDS_genes, bindingRegions = c(-within_genes_kb, within_genes_kb),
    featureType = "Exon",
    PeakLocForDistance = "end",
    FeatureLocForDistance = "start",
    output = "nearestLocation")  
  

  stacked = c(overlapping_ap,upstream_of_peak_togene_start_ap,upstream_of_peak_togene_end_ap,downstream_of_peak_togene_start_ap,downstream_of_peak_togene_end_ap)
  stacked = stacked[!is.na(stacked$insideFeature)]
  stacked = stacked[ order(stacked$peak,stacked$shortestDistance)]
  stacked_nr = stacked[ unique_indexes(stacked$peak)]
  stacked_nr$insideFeature = as.character(stacked_nr$insideFeature)
  stacked_nr[stacked_nr$fromOverlappingOrNearest != 'Overlapping' & stacked_nr$shortestDistance > 5000]$insideFeature <- 'unmapped'
  stacked_nr[stacked_nr$insideFeature == 'unmapped']$feature <- NA
  ap = stacked_nr
  
  # The output is X00001.WBID. The peak ID (in the name column) is more useful.
  names(ap) <- ap$name
  # these fields got converted to character during bigbed output
  # TODO- just do that conversion during the export
  ap$k4weights = as.numeric(ap$k4weights)
  ap$k11weights = as.numeric(ap$k11weights)
  ap$LE_nonNormed = as.numeric(ap$LE_nonNormed)
  ap$L1_nonNormed = as.numeric(ap$L1_nonNormed)
  ap$L3_nonNormed = as.numeric(ap$L3_nonNormed)
  ap$L3_std = as.numeric(ap$L3_std)
  ap$L1_std = as.numeric(ap$L1_std)
  ap$LE_std = as.numeric(ap$LE_std)
  #ap$peak = as.integer(ap$peak) # this are the string "ELT2peak..."
  
  # unmapped +/- 5Kb
  nappy = ap[is.na(ap$feature)]
  # mapped
  ap = ap[!is.na(ap$feature)]
  # reduce list down to the minumum absolute value of distancetoFeature
  ap.ordered = ap[ order(ap$name, abs(ap$distancetoFeature) )]
  ap.unique = ap[ match(unique(ap.ordered$name), ap.ordered$name)]
  ap.rejoined = c(ap.unique, nappy)
  ap.rejoined = ap.rejoined[ order(ap.rejoined$name)]
  
  ap = ap.rejoined
  insideFeatureLabels = as.character(ap$insideFeature)
  insideFeatureLabels[is.na(ap$feature)] <- sprintf("unmapped ± %dKb",within_genes_kb) 
  insideFeatureLabels[ap$shortestDistance > within_genes_kb * 1000] <- sprintf("unmapped ± %dKb",5) 
  ap$insideFeature <- as.factor(insideFeatureLabels)
  # does not differ between clusters, with overlapStart and upstream accounting for 60-70% of the annotation types (followed by 'inside', then 'downstream')
  round(table(ap$k4cluster, ap$insideFeature)/apply(table(ap$k4cluster, ap$insideFeature), 1, sum),3)*100
  
  # ap broken down by kclust number
  ap_0 = ap[ap$k4cluster ==0]
  ap_1 = ap[ap$k4cluster ==1]
  ap_2 = ap[ap$k4cluster ==2]
  ap_3 = ap[ap$k4cluster ==3]
  ap_4 = ap[ap$k4cluster ==4]
  saveRDS(ap, "annotatedPeaks.rds")
  return(
    list(ap=ap, 
         stacked_nr=stacked_nr,
         bycluster=list(
           ap_0=ap_0,
           ap_1=ap_1,
           ap_2=ap_2,
           ap_3=ap_3,
           ap_4=ap_4),
         all_CDS_genes=all_CDS_genes
         )
    )
}