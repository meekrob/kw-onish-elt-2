# Here is my character vector of IDs. Needs to be included in the 'attributes' argument of biomaRt::getBM
length(argument_ids)

# dependancies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("ensembldb", quietly = TRUE)) {
    BiocManager::install("ensembldb")
}
library(ensembldb)

if (!requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt")
}
library(biomaRt)
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
, annot = annFUN.gene2GO # YES!
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


###
### subclasses versus the main class
###
length(ap_0)
length(ap_1)
length(ap_2)
length(ap_3)
length(ap_4)

ap_0_unique_wbid = unique(ap_0$feature)
ap_1_unique_wbid = unique(ap_1$feature)
ap_2_unique_wbid = unique(ap_2$feature)
ap_3_unique_wbid = unique(ap_3$feature)
ap_4_unique_wbid = unique(ap_4$feature)

# now the master set is the previous argument set: only the genes mapped to peaks instead of the whole genome
allpeaks.genes = argument_ids
k0_geneList = factor(as.integer(allpeaks.genes %in% ap_0_unique_wbid))
names(k0_geneList) <- allpeaks.genes

BP_0.go = new("topGOdata", ontology='BP'
, allGenes = k0_geneList
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)

k0_annotated_BP = GOSummary(BP_0.go)
k0_annotated_BP$elim = as.numeric(k0_annotated_BP$elim)

View(k0_annotated_BP)

k1_geneList = factor(as.integer(allpeaks.genes %in% ap_1_unique_wbid))
names(k1_geneList) <- allpeaks.genes

BP_1.go = new("topGOdata", ontology='BP'
, allGenes = k1_geneList
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)

k1_annotated_BP = GOSummary(BP_1.go)
k1_annotated_BP$elim = as.numeric(k1_annotated_BP$elim)

View(k1_annotated_BP)

k2_geneList = factor(as.integer(allpeaks.genes %in% ap_2_unique_wbid))
names(k2_geneList) <- allpeaks.genes

BP_2.go = new("topGOdata", ontology='BP'
, allGenes = k2_geneList
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)

k2_annotated_BP = GOSummary(BP_2.go)
k2_annotated_BP$elim = as.numeric(k2_annotated_BP$elim)

View(k2_annotated_BP)


k3_geneList = factor(as.integer(allpeaks.genes %in% ap_3_unique_wbid))
names(k3_geneList) <- allpeaks.genes

BP_3.go = new("topGOdata", ontology='BP'
, allGenes = k3_geneList
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)

k3_annotated_BP = GOSummary(BP_3.go)
k3_annotated_BP$elim = as.numeric(k3_annotated_BP$elim)

View(k3_annotated_BP)

k4_geneList = factor(as.integer(allpeaks.genes %in% ap_4_unique_wbid))
names(k4_geneList) <- allpeaks.genes

BP_4.go = new("topGOdata", ontology='BP'
, allGenes = k4_geneList
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)

k4_annotated_BP = GOSummary(BP_4.go)
k4_annotated_BP$elim = as.numeric(k4_annotated_BP$elim)

View(k4_annotated_BP)
