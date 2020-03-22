if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
library(biomaRt)

paramart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)
all_transcripts = getBM(mart = paramart, 
                     filter=c("species_id_1010","chromosome_name"), 
                     value=list(species_id_1010="caelegprjna13758",chromosome_name="I,II,III,IV,V,X"), #36 Mt entries ignored
                     attributes= c("chromosome_name", 
                                   "start_position",
                                   "end_position",
                                   "wbps_transcript_id",
                                   "strand",
                                   # score would go here in BED
                                   "genomic_coding_start",
                                   "genomic_coding_end",
                                   # fields to get the blockCount, blockSizes and blockStarts from
                                   "exon_chrom_start",
                                   "exon_chrom_end",
                                   "rank",
                                   # extra annotations to retain
                                   "wbps_gene_id",
                                   "external_gene_id",
                                   "description",
                                   "entrezgene_id",
                                   "entrezgene_name",
                                   "gene_biotype",
                                   "wormbase_locus",
                                   "wormbase_gseq",
                                   "refseq_peptide",
                                   "refseq_mrna"
                                   )
              )
all_transcripts$strand = ifelse(all_transcripts$strand==1,'+','-')
all_transcripts$blockSize = all_transcripts$exon_chrom_end - all_transcripts$exon_chrom_start
all_transcripts$blockStart = all_transcripts$exon_chrom_start - all_transcripts$start_position
all_transcripts$chromosome_name = paste("chr", all_transcripts$chromosome_name, sep='') 
table(all_transcripts$chromosome_name)

all_transcripts$exon_chrom_start = NULL
all_transcripts$exon_chrom_end = NULL
all_transcripts[na_CDS,'genomic_coding_start'] = all_transcripts[na_CDS, 'start_position']
all_transcripts[na_CDS,'genomic_coding_end'] = all_transcripts[na_CDS, 'start_position']
all_transcripts[all_transcripts == ''] <- NA
# we'll use the wormbase_locus ID as the name, unless it's unset, then use the transcript ID
all_transcripts$wormbase_locus[ is.na(all_transcripts$wormbase_locus)] <- all_transcripts[is.na(all_transcripts$wormbase_locus),'wbps_transcript_id']

# colors for types of genes
gene_colors = t(col2rgb(RColorBrewer::brewer.pal(9, 'Paired')))
gene_colors = rbind(c(0,0,0), c(128,128,128), gene_colors) # add protein and pseudogene
biotypes = unique(all_transcripts$gene_biotype)
# paste the output of this line into 
# it will produce an extra comma in the front. Delete it.
cat("colors = {", 
    sprintf(",\n'%s': '%d,%d,%d'", 
            rownames(gene_colors), 
            gene_colors[,1],
            gene_colors[,2],
            gene_colors[,3]),"}",sep='')
# colors = {
#   'protein_coding': '0,0,0',
#   'pseudogene': '128,128,128',
#   'tRNA': '166,206,227',
#   'miRNA': '31,120,180',
#   'ncRNA': '178,223,138',
#   'rRNA': '51,160,44',
#   'snRNA': '251,154,153',
#   'snoRNA': '227,26,28',
#   'lincRNA': '253,191,111',
#   'piRNA': '255,127,0',
#   'antisense_RNA': '202,178,214'}

# paste the output of this command into the html detail page
cat("<table>", 
    sprintf("<tr><td style='background-color:rgb(%d,%d,%d)'>%s</td><td>%s</td></tr>\n", 
            gene_colors[,1],
            gene_colors[,2],
            gene_colors[,3],
            rownames(gene_colors),
            rownames(gene_colors)),
            "</table>",
            sep='')

# <table><tr><td style='background-color:rgb(0,0,0)'>protein_coding</td><td>protein_coding</td></tr>
#   <tr><td style='background-color:rgb(128,128,128)'>pseudogene</td><td>pseudogene</td></tr>
#   <tr><td style='background-color:rgb(166,206,227)'>tRNA</td><td>tRNA</td></tr>
#   <tr><td style='background-color:rgb(31,120,180)'>miRNA</td><td>miRNA</td></tr>
#   <tr><td style='background-color:rgb(178,223,138)'>ncRNA</td><td>ncRNA</td></tr>
#   <tr><td style='background-color:rgb(51,160,44)'>rRNA</td><td>rRNA</td></tr>
#   <tr><td style='background-color:rgb(251,154,153)'>snRNA</td><td>snRNA</td></tr>
#   <tr><td style='background-color:rgb(227,26,28)'>snoRNA</td><td>snoRNA</td></tr>
#   <tr><td style='background-color:rgb(253,191,111)'>lincRNA</td><td>lincRNA</td></tr>
#   <tr><td style='background-color:rgb(255,127,0)'>piRNA</td><td>piRNA</td></tr>
#   <tr><td style='background-color:rgb(202,178,214)'>antisense_RNA</td><td>antisense_RNA</td></tr>
#   </table>

write.table(all_transcripts[all_transcripts$blockSize > 0,], "all_exons.table", sep="\t", col.names=T, row.names=F,quote=F)
system("python3 collapse_exons.py all_exons.table wbps_transcript_id> everything.collapsed")
system("python3 process_collapsed.py everything.collapsed > everything.bedPlus")
system("/Users/david/bin/UCSC_userApps/bedSort everything.bedPlus everything.bedPlus")
system("/Users/david/bin/UCSC_userApps/bedToBigBed everything.bedPlus chrom.sizes WS271.bb -type=bed12+10 -tab -as=WS271plus.as")
