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

listEnsembl() # ensembl      Ensembl Genes 96
ensembl <- useEnsembl(biomart = "ensembl")
listDatasets(mart = ensembl) # Use number 29: celegans_gene_ensembl
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset = "celegans_gene_ensembl")
# Get TSS
ce_TSS = getAnnotation(mart, featureType = "TSS")

source('scripts/tenthousandPeaksAsGRanges.R') # defines "peaks"
annotated_peaks = annotatePeakInBatch(peaks,mart, featureType = "TSS", AnnotationData = ce_TSS)
wbid=unlist(strsplit(names(annotated_peaks), "[.]"))[seq(2,2*length(annotated_peaks),by=2)]
annotated_peaks$WBID = wbid
