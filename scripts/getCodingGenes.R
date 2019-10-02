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
