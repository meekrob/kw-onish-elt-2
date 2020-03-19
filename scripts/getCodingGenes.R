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

PKG="ChIPpeakAnno"
if (!requireNamespace(PKG, quietly = TRUE)) { BiocManager::install(PKG) }
library(PKG, character.only=T)
getCodingGenes = function(peaks){
  # Get only the protein-coding genes. Ranges are comprehensive across splice variants.
  paramart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)
  genes_coding = getBM(mart = paramart, 
                       filter=c("species_id_1010", 
                                "biotype"), 
                       value=list(species_id_1010="caelegprjna13758", 
                                  biotype="protein_coding"), 
                       attributes = c('wbps_gene_id',
                                      'chromosome_name', 
                                      'start_position', 
                                      'end_position', 
                                      'strand',
                                      'wormbase_gseq'))
  # exclude mitochondria
  genes_coding_noMT = genes_coding[genes_coding$chromosome_name != 'MtDNA',]
  
  bedformat = genes_coding_noMT[,c('chromosome_name','start_position','end_position','strand','wbps_gene_id')]
  bedformat$strand = with(bedformat, ifelse(strand==1, '+','-'))
  bedformat = bedformat[with(bedformat, order(chromosome_name, start_position)), ] # WS271 as of this writing
  write.table(bedformat, "celegans_genes.WS271.bed",col.names = F,row.names = F,quote=F, sep="\t")
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
  ap=annotatePeakInBatch(peaks, AnnotationData=all_CDS_genes)
  
  # a pie chart with the breakdown of how the annotation happened
  pie_labels = paste0(names(table(ap$insideFeature)), rep(" (",5), table(ap$insideFeature), rep(")",5))
  pie(table(ap$insideFeature), labels=pie_labels, main="ChIPpeakAnno::annotatePeakInBatch: \"nearestLocation\"", sub="paramart: wbps_gene_id, biotype= protein_coding")
  # does not differ between clusters, with overlapStart and upstream accounting for 60-70% of the annotation types (followed by 'inside', then 'downstream')
  round(table(ap$k4cluster, ap$insideFeature)/apply(table(ap$k4cluster, ap$insideFeature), 1, sum),3)*100
  
  
  ap.wbid=unlist(strsplit(names(ap), "[.]"))[seq(2,2*length(ap),by=2)]
  unique.ap.wbid = unique(ap.wbid)
  
  # ap broken down by kclust number
  ap_0 = ap[ap$k4cluster ==0]
  ap_1 = ap[ap$k4cluster ==1]
  ap_2 = ap[ap$k4cluster ==2]
  ap_3 = ap[ap$k4cluster ==3]
  ap_4 = ap[ap$k4cluster ==4]
  
  return(
    list(ap=ap, 
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