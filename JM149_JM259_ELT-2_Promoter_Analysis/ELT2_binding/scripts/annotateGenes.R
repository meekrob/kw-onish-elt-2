# set up environment, install packages if needed
PKG="BiocManager" # get the Bioc installer to get the rest
if (!requireNamespace(PKG, quietly = TRUE)) { install.packages(PKG) }
library(PKG, character.only=T)

PKG="ChIPpeakAnno"
if (!requireNamespace(PKG, quietly = TRUE)) { BiocManager::install(PKG) }
library(PKG, character.only=T)

PKG="GenomicFeatures"
if (!requireNamespace(PKG, quietly = TRUE)) { BiocManager::install(PKG) }
library(PKG, character.only=T)

PKG="biomaRt"
if (!requireNamespace(PKG, quietly = TRUE)) { BiocManager::install(PKG) }
library(PKG, character.only=T)

PKG="BSgenome"
if (!requireNamespace(PKG, quietly = TRUE)) { BiocManager::install(PKG) }
library(PKG, character.only=T)

# connect to a databases
paramart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)
mart = useMart(biomart = "ensembl", dataset = "celegans_gene_ensembl")

# this returns transcript accessions with a suffix added, which will require processing to get the Wormbase ID
ExonPlusUtr.ce11 = getAnnotation(mart=mart,featureType="ExonPlusUtr")

# these can be mapped to Wormbase IDs, see "addGeneIDs" below
TSS.ce11 = getAnnotation(mart=mart,featureType="TSS")

# make "annotated" by running ChIPpeakAnno::annotatePeakInBatch on a set of GRanges
stopifnot("LE_IDR" %in% ls())

# still looking for the best way to do this...
annotatePeakInBatch(LE_IDR, paramart, featureType="TSS",AnnotationData = TSS.ce11, output="both", maxgap=5000) 
-> annotated
annotatePeakInBatch(LE_IDR, featureType="TSS",AnnotationData = ExonPlusUtr.ce11, output="both", maxgap=5000) 
-> annotated

# add Wormbase ID mappings
addGeneIDs(annotated, mart=paramart, feature_id_type="wbps_transcript_id", IDs2Add=c("wbps_gene_id"))->annotated
