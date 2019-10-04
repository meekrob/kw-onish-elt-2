# Here is my character vector of IDs. Needs to be included in the 'attributes' argument of biomaRt::getBM
length(argument_ids)

# dependancies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("topGO", quietly = TRUE)) {
    BiocManager::install("topGO")
}
library(topGO)


# connect to public database for annotations and geneset
paramart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

#   Get the specs from https://parasite.wormbase.org/biomart/martview/ 
#   You will need to add transcript IDs to your output columns (under attributes) if that's 
# what you're starting with.
WORMGO=biomaRt::getBM(mart = paramart, filter="species_id_1010", value="caelegprjna13758",
   attributes = c( "wbps_gene_id", 
        "external_gene_id", 
        "go_accession", 
        "go_name_1006", 
        "go_linkage_type"))


# Remove unannotated. I don't know how this affects the enrichment, it was in an example.
unannotated_genes = WORMGO[WORMGO$go_accession == '',]
WORMGO = WORMGO[WORMGO$go_accession != '',]

# create an object where you can access all the GO terms that are assigned to a specific gene
geneID2GO <- by(WORMGO$go_accession, WORMGO$wbps_gene_id, function(x) as.character(x))

# Master set of genes (only the ones that have annotation).
# nonredundant list of genes in genome
all.genes <- unique(as.character(WORMGO$wbps_gene_id)) # could also get it via: names(geneID2GO)

# a boolean table, one entry for each of the genes in the master set. 1 if its in our argument set, 0 if not.
geneList = factor(as.integer(all.genes %in% argument_ids))
names(geneList) = all.genes

# A helper function to do all the statistical testing and put them together. I have been using 'elim', but comparing the 
# result to other statistics.
GOSummary<- function(GOdata) {
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFis <- getSigGroups(GOdata, test.stat)
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(GOdata, test.stat)
  test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test",cutOff = 0.01)
  resultElim <- as.numeric( getSigGroups(GOdata, test.stat) ) # Field of interest for me. as.numeric will put in NAs for strings like '< 1e-30' though.
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(GOdata, test.stat)
  l <- list(classic = resultFis, KS = resultKS, elim = resultElim,weight = resultWeight)
  return(GenTable(object=GOdata, weight=l$weight, classic=l$classic, elim=l$elim, KS=l$KS, orderBy="weight",ranksOf = "classic", topNodes = 50))
}

BP.go = new("topGOdata", ontology='BP'
, allGenes = geneList
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)

MF.go = new("topGOdata", ontology='MF'
, allGenes = geneList
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)

CC.go = new("topGOdata", ontology='CC'
, allGenes = geneList
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)

# cuts things down to top 50 terms. See 'topNodes' argument in 'GenTable' called by 'GOSummary'.
annotated_BP = GOSummary(BP.go)
annotated_MF = GOSummary(MF.go)
annotated_CC = GOSummary(CC.go)

# now I look at the results like
View(annotated_BP)
# and maybe filter via a threshold
annotated_BP_sig = annotated_BP[annotated_BP$elim <= .1]
