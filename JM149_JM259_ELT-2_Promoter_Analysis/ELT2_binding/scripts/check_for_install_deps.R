# standard packages
pkgs = c(
  "BiocManager",
  "reshape2",
  "stringr",
  "pheatmap",
  "ggplot2",
  "ggExtra",
  "ggrepel",
  "RColorBrewer",
  "VennDiagram",
  "dplyr",
  "knitr"
)
need.pkgs = pkgs[! (pkgs %in% installed.packages()[,"Package"])]
if (length(need.pkgs)) {install.packages(need.pkgs, update=F) }

# Bioconductor packages
pkgs = c("biomaRt",
         "GenomicRanges", 
         "ChIPpeakAnno",
         "topGO", 
         "plyranges", # tidy the GenomicRanges objects
         "UpSetR")

need.pkgs = pkgs[! (pkgs %in% installed.packages()[,"Package"])]
if (length(need.pkgs)) {
  library(BiocManager)
  BiocManager::install(need.pkgs,update=F)
}
