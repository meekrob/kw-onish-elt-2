## GetWSpromoters
# limitations:
# several transcripts differ in bioMart by a small number of positions,
# leading to potentially crazy overrepresentation by them. Example:

# chrI    10232   11232   -       WBGene00022277  Y74C9A.3        homt-1
# chrI    10495   11495   +       WBGene00022276  Y74C9A.2        nlp-40
# chrI    10524   11524   +       WBGene00022276  Y74C9A.2        nlp-40
# chrI    10606   11606   -       WBGene00022277  Y74C9A.3        homt-1
# chrI    14103   15103   +       WBGene00022276  Y74C9A.2        nlp-40
# chrI    21127   22127   -       WBGene00022278  Y74C9A.4        rcor-1
# chrI    26643   27643   -       WBGene00022278  Y74C9A.4        rcor-1
# chrI    26781   27781   -       WBGene00022278  Y74C9A.4        rcor-1
# chrI    26782   27782   -       WBGene00022278  Y74C9A.4        rcor-1


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

getPromoters = function(upstream = 1000)
{
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
      'chromosome_name',
      'transcript_start',
      'transcript_end',
  #    'start_position',
  #    'end_position',
      'strand',
      'wbps_gene_id',
      'wormbase_gseq',
      'wormbase_locus'
    )
  )
  genes_coding = genes_coding %>% filter(chromosome_name != "MtDNA")
  genes_coding = genes_coding %>% mutate(chromosome_name = sprintf("chr%s", chromosome_name))
  genes_coding$promoter_start = ifelse(genes_coding$strand == 1, 
                                         genes_coding$transcript_start - upstream, 
                                         genes_coding$transcript_end)
  genes_coding$promoter_end = ifelse(genes_coding$strand == 1, 
                                         genes_coding$transcript_start , 
                                         genes_coding$transcript_end + upstream)
  promoters = with(genes_coding, 
                   data.frame(chrom  = chromosome_name,
                              start  = promoter_start, 
                              end    = promoter_end,
                              strand = ifelse(strand==1, "+", "-"),
                              wbid = wbps_gene_id,
                              wormbase_gseq = wormbase_gseq,
                              wormbase_locus = wormbase_locus))
  promoters
}
# usage
#write.table(getPromoters(1000),"promotersWS1Kb.bed", quote=F, sep="\t", row.names=F,col.names=F)
# in bash, you must :
# sort -k1,1 -k2,2n promotersWS1Kb.bed | uniq > tmp && mv tmp promotersWS1Kb.bed