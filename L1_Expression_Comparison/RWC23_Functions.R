myLibs <- function(){
  library(biomaRt)
  library(tidyverse)
  library(readxl)
  library(ComplexHeatmap)
  library(matrixStats)
  library(pheatmap)
  library(RVAideMemoire)
  library(dendextend)
  library(binom)
  library(circlize)
}



# matrix_select() subsets a matrix based on vactor of gene names.
# Pleasure ensure gene naming is consistent between files.
# It is best to use WBGeneIDs for this process

matrix_select <- function(count_matrix, gene_subset_vector){
  nameselect <- rownames(count_matrix) %in% gene_subset_vector
  nameMatrix <- count_matrix[nameselect,]
  return(nameMatrix)
}

# row_normalize_matrix() performs a row normalization by subtracting the row mean from each value in the row.

row_normalize_matrix <- function(count_matrix){
  namevarRowNormalized <- count_matrix - rowMeans(count_matrix)
  return(namevarRowNormalized)
}

# row_normalize_matrix_cutoff() subsets the count matrix based on value in variance_cutoff
# It then performs a row normalization by subtracting the row mean from each value in the row.

row_normalize_matrix_cutoff <- function(count_matrix, variance_cutoff){
  varselect <- rowVars(count_matrix) >= variance_cutoff
  count_matrix <- count_matrix[varselect, ]
  namevarRowNormalized <- count_matrix - rowMeans(count_matrix)
  return(namevarRowNormalized)
}

# row_zscore_matrix() assigns a row Z score to each row value. 

row_zscore_matrix <- function(count_matrix){
  namevarZscore <- (count_matrix - rowMeans(count_matrix))/rowSds(count_matrix)
  return(namevarZscore)
  
}

# row_zscore_matrix_cutoff() assigns a row Z score to each row value. 
# Then uses variance to subset data.

row_zscore_matrix_cutoff <- function(count_matrix, variance_cutoff){
  varselect <- rowVars(count_matrix) >= variance_cutoff
  count_matrix <- count_matrix[varselect, ]
  namevarZscore <- (count_matrix - rowMeans(count_matrix))/rowSds(count_matrix)
  return(namevarZscore)
  
}

# id2name() replaces WBGeneIDs in The input count matrix to 
# gene names for easier interpretation of genes on a heatmap.
# This function requires access to the internet for proper funcitoning.

id2name <- function(count_matrix){
  paramart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)
  id2name.df <- suppressMessages(biomaRt::getBM(mart = paramart, 
                                                filter=c("wbps_gene_id"), 
                                                value=rownames(count_matrix), 
                                                attributes = c('wbps_gene_id','external_gene_id')
  ))
  rownames(count_matrix) <- id2name.df$external_gene_id[ 
    match(rownames(count_matrix), id2name.df$wbps_gene_id)
    ]
  return(count_matrix)
}

# myPheatmap() is a custom pheatmap plotting function

myPheatmap <- function(count_matrix, title, rowspace){
  pheatmap(count_matrix, 
           cluster_cols = FALSE, 
           cluster_rows = TRUE, 
           show_rownames = FALSE, 
           border_color = NA,
           main = title,
           cutree_rows = rowspace)
}

mysmallPheatmap <- function(count_matrix, title, rowspace){
  pheatmap(count_matrix, 
           cluster_cols = FALSE, 
           cluster_rows = TRUE, 
           show_rownames = TRUE, 
           border_color = NA,
           main = title,
           cellheight = 10,
           cutree_rows = rowspace)
}



# Use the binomial test to calculate a pvalue from a contingency table

ctable_binom <- function(ctable, alt, p = 0.05){
  df <- data.frame()
  x <- 1
  for (i in rownames(ctable)){
    xbinom = as.numeric(ctable[i,1]) # number bound in a set
    nbinom = as.numeric(ctable[i,1] + ctable[i,2]) # total genes in a set
    pval <- binom.test(xbinom,nbinom,proportion, alternative = alt)$p.value # calculate the pvalue
    conf95 <- binom.test(xbinom,nbinom,proportion, alternative = alt)$conf.int
    toappend <- data.frame(Set = i, pval = pval, conf.lower = conf95[1], conf.upper = conf95[2], stringsAsFactors = FALSE)
    df <- bind_rows(df, toappend)
    x <- x + 1
  }
  print(df %>% mutate(bool = pval < p))
}

# RNA_heatmap() uses ComplexHeatmap to make a heatmap 

RNA_heatmap <- function(mat, split = NULL){
  Heatmap(
    mat,
    name = "elt2D-elt7D\nRNAseq",
    col = colorRampPalette(c("cyan", "black", "yellow"))(1000),
    cluster_columns = FALSE,
    column_split = factor(
      c(
        rep("wt", 4),
        rep("elt7D", 3),
        rep("elt2D", 4),
        rep("elt7Delt2D", 3)
      ),
      levels = c("wt", "elt7D", "elt2D", "elt7Delt2D")
    ),
    split = split,
    clustering_distance_rows = "spearman",
    clustering_method_rows = "complete",
    show_row_names = FALSE,
    show_column_names = TRUE,
    row_names_gp = gpar(cex = 0.2),
    column_names_gp = gpar(cex = 0.4),
    heatmap_legend_param = list(color_bar = "continuous"),
    bottom_annotation = HeatmapAnnotation(foo = anno_block(
      labels = c("wt", "elt7D", "elt2D", "elt7Delt2D"),
      gp = gpar(border = NA, lty = "blank"),
    ))
  )
}

# RNA_heatmap2() is a variation with better plotting aestetics

RNA_heatmap2 <- function(mat, column_split = NULL, row_title = NULL, row_split = NULL, split = NULL, right_annotation = NULL){
  Heatmap(
    mat,
    name = "elt2D-elt7D\nRNAseq",
    col = colorRampPalette(c("cyan", "black", "yellow"))(1000),
    cluster_columns = FALSE,
    clustering_distance_rows = "spearman",
    clustering_method_rows = "complete",
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_labels = column_labels[colnames(mat)],
    column_names_gp = gpar(cex = 0.7),
    heatmap_legend_param = list(color_bar = "continuous"),
    show_row_dend = FALSE,
    row_title_rot = 0,
    row_split = row_split,
    cluster_row_slices = FALSE,
    # row_title = row_title,
    column_title = NULL,
    column_split = column_split,
    right_annotation = right_annotation,
    bottom_annotation = HeatmapAnnotation(
      foo = anno_block(
        labels = c("WT", "elt-7(tm840)", "elt-2(ca15)", "elt-7(tm840);\nelt-2(ca15)"),
        labels_gp = gpar(cex = .8),
        gp = gpar(border = NA, lty = "blank")
      ),
      foo2 = anno_block(gp = gpar(fill = "black"), height = unit(0.5, "mm"))
    ),
    # left_annotation = rowAnnotation(foo = anno_block(
    #   labels = c("SET1", "SET2", "SET3", "SET4", "SET5", "SET6"),
    #   labels_rot = 0,
    #   gp = gpar(border = NA, lty = "blank", cex = 0.4)
    # ))
  )
}


# elt2_l1_row_annotation() should be used in conjuntion with RNA_heatmap().
# It will add row annotations to RNA_heatmap for ELT2 binding.
# The `df` must be the same size and in the same row order as `mat` for RNA_heatmap().

elt2_l1_row_annotation <- function(df){
  rowAnnotation(L1bound = df$elt2_detected_in_L1,
                col = list(L1bound = c(
                  "Direct" = "green", "Indirect" = "black"
                )))
}

# binding_cluster_row_annotation() will add row annotations based on what binding cluster the gene is assocaited with.
# The `df` must be the same size and in the same row order as `mat` for RNA_heatmap().

binding_cluster_row_annotation <- function(df){
  rowAnnotation(
    Embryo_Specific = df$Embryo_Specific,
    col = list(Embryo_Specific = c(
      "absent" = "white", "present" = "#7570B3"
    )),
    border = TRUE
  ) +
    rowAnnotation(Larval = df$Larval,
                  col = list(Larval = c(
                    "absent" = "white", "present" = "#1B9E77"
                  )),
                  border = TRUE) +
    rowAnnotation(
      Increasing = df$Increasing,
      col = list(Increasing = c(
        "absent" = "white", "present" = "#E7298A"
      )),
      border = TRUE
    ) +
    rowAnnotation(
      L3_High = df$L3_High,
      col = list(L3_High = c(
        "absent" = "white", "present" = "#D95F02"
      )),
      border = TRUE
    ) +
    rowAnnotation(
      Not_Changing = df$Not_Changing,
      col = list(Not_Changing = c(
        "absent" = "white", "present" = "#505050"
      )),
      border = TRUE
    # ) +
    # rowAnnotation(
    #   Not_Bound = df$Not_Bound,
    #   col = list(Not_Bound = c(
    #     "absent" = "white", "present" = "black"
    #   )),
    #   border = TRUE
    )
}

# make_cluster_annotation() logs the number of binding sites per gene (row) that has an ELT2 binding pattern

make_cluster_annotation <- function(input_matrix, binding_matrix){
  df <- input_matrix %>%
    as.data.frame.matrix() %>%
    rownames_to_column() %>%
    left_join(rownames_to_column(binding_matrix), by = "rowname") %>%
    dplyr::select(rowname, all_of(elt2_cluster_names), any_of("Not_Bound")) %>%
    replace_na(list(
      Embryo_Specific = 0,
      Larval = 0,
      Increasing = 0,
      L3_High = 0,
      Not_Changing = 0
    )) %>% 
    dplyr::rename(WBGeneID = "rowname")
  df
}

# make_cluster_binary_annotation() will convert number of binding sites per ELT2 binding pattern
# to binary present/absent list. Used in conjunction with make_cluster_annotation()
# `input_matrix` is the output of make_cluster_annotation()

make_cluster_binary_annotation <- function(input_matrix){
  input_matrix %>%
    mutate(
      Embryo_Specific = if_else(
        condition = input_matrix$Embryo_Specific == 0,
        true = "absent",
        false = "present"
      )
    ) %>%
    mutate(Larval = if_else(
      condition = input_matrix$Larval == 0,
      true = "absent",
      false = "present"
    )) %>%
    mutate(
      Increasing = if_else(
        condition = input_matrix$Increasing == 0,
        true = "absent",
        false = "present"
      )
    ) %>%
    mutate(
      L3_High = if_else(
        condition = input_matrix$L3_High == 0,
        true = "absent",
        false = "present"
      )
    ) %>%
    mutate(
      Not_Changing = if_else(
        condition = input_matrix$Not_Changing == 0,
        true = "absent",
        false = "present"
      )
    # ) %>%
    # mutate(
    #   Not_Bound = if_else(
    #     condition = input_matrix$Not_Bound == 0,
    #     true = "absent",
    #     false = "present"
    #   )
    )
  }

row_scale <- function(mat){
  scale_mat <-t(apply(unlist(mat), 1, scale))
  colnames(scale_mat) <- colnames(mat)
  rownames(scale_mat) <- rownames(mat)
  scale_mat
}

# GOI_annotate_heatmap() will identify indicies for genes of interest (GOI) from a gene expression dataset
# the gene expression dataset must be a dataframe with a column that contains gene names or WBGeneIDs
# Input a vector of gene names and corresponding data frame column vector (using "$" operator)

GOI_annotate_heatmap <- function(gene_names, rna_df_column){
  index_list <- sapply(gene_names, function(x){
    which(rna_df_column == x)
  }
  )
  data.frame(index = index_list) %>% rownames_to_column(var = "name")
}

# myPDFplot(), wrapper to save ComplexHeatmap images to a PDF file

myPDFplot <- function(plot, name, height, width, plotdir = plotdir) {
  pdf(
    paste(plotdir,
          name,
          "_",
          today(),
          ".pdf",
          sep = ""),
    height = height,
    width = width
  )
  print(plot)
  dev.off()
}

myggsave <- function(plot, name, height = NA, width = NA, plotdir = plotdir) {
    ggsave(
      filename = paste(plotdir,
            name,
            "_",
            today(),
            ".pdf",
            sep = ""),
      plot = plot,
      height = height,
      width = width
      )
}
