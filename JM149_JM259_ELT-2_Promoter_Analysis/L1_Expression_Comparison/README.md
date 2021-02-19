# RWC23_elt2_regulated_genes

The goal of this analysis is to better understand developmental ChIP-seq data of the C. elegans intestine master regulator ELT-2.

There are three parts to this analysis

1. Integrate ELT-2 ChIP-seq with ELT-2/EL-7 differential expression
2. Gather publically available C. elegans intestine datasets
3. Integrate ELT-2 ChIP-seq with Time-resolved C. elegans gene expression

## Directory structure

- Each of the analysis parts listed above are contained within a directory
- Corresponding input, output and plots are contained within a respective directory
- Analysis is performed using `.Rmd` Rmarkdown files
- Analysis with plotting components have a companion knitr RMarkdown PDF file
- Custom functions are stored in the file `RWC23_Functions.R`


# Data sources

## Differential expression data

This data was originally published in the following paper:

> Dineen A, Osborne Nishimura E, Goszczynski B, Rothman JH, McGhee JD. Quantitating transcription factor redundancy: The relative roles of the ELT-2 and ELT-7 GATA factors in the C. elegans endoderm. Dev Biol. January 2018. doi:10.1016/j.ydbio.2017.12.023

The file `Table_S1_Raw_Read_Counts.xlsx` was downloaded and converted to a `.csv` file.

## Public C. elegans Intestine Expression

> Tintori SC, Nishimura EO, Golden P, Lieb JD, Goldstein B. A Transcriptional Lineage of the Early C . elegans Embryo. Dev Cell. 2016;38(4):430-444. doi:10.1101/047746

> Hashimshony T, Feder M, Levin M, Hall BK, Yanai I. Spatiotemporal transcriptomics reveals the evolutionary history of the endoderm germ layer. Nature. 2015;519(7542):219-222. doi:10.1038/nature13996

> Pauli F, Liu Y, Kim Y a, Chen P-J, Kim SK. Chromosomal clustering and GATA transcriptional regulation of intestine-expressed genes in C. elegans. Development. 2006;133(2):287-295. doi:10.1242/dev.02185

> Spencer, W. C., Zeller, G., Watson, J. D., Henz, S. R., Watkins, K. L., McWhirter, R. D., … Miller, D. M. A spatial and temporal map of C. elegans gene expression. Genome Research, 21(2), 325–41. https://doi.org/10.1101/gr.114595.110


## ChIP-seq

This data was originally published in the following paper:

> Kudron MM, Victorsen A, Gevirtzman L, et al. The modern resource: genome-wide binding profiles for hundreds of Drosophila and Caenorhabditis elegans transcription factors. Genetics. 2018;208(3):937-949. doi:10.1534/genetics.117.300657

ELT-2 ChIP-seq raw reads were processed by David King, and generated the file `L1_1_L1_2.IDR_0.05.narrowPeak`. See the following reposity for more information: https://github.com/meekrob/onish_ChIP_R_Analysis

## Issues to note

- RNA-seq data aligned with ce10 assembly with WS235 annotation
- Chip-seq data aligned to ce11 assembly
