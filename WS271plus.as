table WS271_plus_xref
"WS271 with multiple other IDs extracted from parasite"
(
string  chrom;		"1: Reference sequence chromosome or scaffold"
uint    chromStart;	"2: Start position of feature on chromosome"
uint    chromEnd;	"3: End position of feature on chromosome"
string  name;		"4: wormbase_locus or Wormbase transcript ID"
uint    score;		"5: Score[unused]"
char[1] strand;		"6: + or - for strand"
uint    thickStart;	"7: Coding region start"
uint    thickEnd;	"8: Coding region end"
uint  	reserved;	"9: Black"
int    blockCount; "10: Number of exons"
int[blockCount]  blockSizes;  "11: Exon sizes"
int[blockCount]  chromStarts; "12: Exon starts"
string wbps_gene_id;     "13: wbps_gene_id"
string external_gene_id; "14: external_gene_id"
string description;      "15: description"
string entrezgene_id;    "16: entrezgene_id"
string entrezgene_name;  "17: entrezgene_name"
string gene_biotype;     "18: gene_biotype"
string wormbase_locus;   "19: wormbase_locus"
string wormbase_gseq;    "20: wormbase_gseq"
string refseq_peptide;   "21: refseq_peptide"
string refseq_mrna;      "22: refseq_mrna"
)
