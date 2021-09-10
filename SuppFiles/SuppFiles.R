library(tidyverse)
library(openxlsx)
# FILE S1
# Intestine development expression cluster assignments
figure3df <- read_csv(file = "./Figure_3_Intestine_Expression_Patterns/02_output/Fig3_Intestine_Development_Cluster_Assignments.csv")
table(figure3df$intestine_dev_cluster)

# ELT-2 ChIP-seq 
elt2_GRange <- readRDS("01_input/annotatedPeaks.rds")
head(elt2_GRange)
elt2_peaks <- as.data.frame(mcols(elt2_GRange))
elt2_peaks <- elt2_peaks %>% dplyr::rename(cluster.description = k4labels, WBGeneID = feature)
elt2_cluster_names <- c("Embryo_Specific",
                        "Larval",
                        "Increasing",
                        "L3_High",
                        "Not_Changing")

elt2_peaks$cluster.description <-
  factor(
    elt2_peaks$cluster.description,
    levels = c(
      "LE-specific",
      "Post-embryonic",
      "Increasing",
      "L3-high",
      "Not-changing or not IDR-passing"
    ),
    labels = elt2_cluster_names
  )

# ELT-2/ELT-7 transcriptional response genes
figure4df <- read_csv(file = "Figure_4_L1_Regulation/02_output/ELT-2_Target_Gene_Class_Assignments.csv")

# merge datasets
elt2_peaks_merged <- elt2_peaks %>% 
  left_join(figure3df, by = "WBGeneID") %>% 
  left_join(figure4df, by = "WBGeneID") %>%
  dplyr::select(peak, WBGeneID, LE_1:L3_std, variance, start_position:fromOverlappingOrNearest, cluster.description, intestine_dev_cluster, target_gene_class = class)

write.xlsx(elt2_peaks_merged, file = "../SuppFiles/FileS1_ELT2_ChIP-seq.xlsx", keepNA = TRUE, na.string = "NA")
# Give better column names

# FILE S2
fileS2 <- read_csv("./Figure_5/02_output/FileS6_ELT2_occupancy_ontology_results.csv") %>% 
  mutate(Ontology.type = case_when(analysis == "go_df" ~ "Gene Ontology", analysis == "phenotype_df" ~ "Phenotype Ontology", analysis == "tissue_df" ~ "Tissue Ontology")) %>%
  dplyr::select(WBGeneID = wbid, Ontology.type, Ontology.name = name, Expected:Q.value, ELT2.occupancy = cluster)

write.xlsx(fileS2, "../SuppFiles/FileS2_ELT2Occupancy_Ontology_Results.xlsx")

# FILE S3

fileS3 <- read_csv("./Figure_3_Intestine_Expression_Patterns/02_output/Fig3_Public_Intestine_RNA_Gene_List.csv")
write.xlsx(fileS3, "../SuppFiles/FileS3_Intestine_Gene_List.xlsx")

# FILE S4

fileS4 <- read_csv("./Figure_3_Intestine_Expression_Patterns/02_output/FileS3_Intestine-specific_time-resolved_data.csv") %>% 
  left_join(figure3df) %>%
  dplyr::select(WBGeneID, Sequence.name, Public_name:spencer_L2, ELT.occupancy = cluster.description,intestine_dev_cluster, emb_510min:YA) %>% 
  filter(!is.na(intestine_dev_cluster)) %>% arrange(WBGeneID)

write.xlsx(fileS4, "../SuppFiles/FileS4_Intestine-specific_time-resolved_RNA-seq.xlsx", overwrite = TRUE)

# FILE S5

fileS5 <- read_csv("./Figure_4_L1_Regulation/02_output/Figure4_WBGeneID_Annotation_Dataframe.csv") %>%
  filter(!(is.na(bound_class))) %>%
  dplyr::select(WBGeneID, Public_name = SYMBOL, TargetGeneClass = bound_class) %>% 
  left_join(read_csv("./Figure_4_L1_Regulation/02_output/Figure4_Dynamic_Counts_BoundOnly.csv"), by = "WBGeneID")
write.xlsx(fileS5, "../SuppFiles/FileS5_TargetGeneClass_RNA-seq.xlsx", keepNA = TRUE, na.string = "NA", overwrite = TRUE)

# FILE S6
fileS6 <- read_csv("./Figure_5/02_output/FileS5_ELT2-Targetgeneclass-Ontology_results.csv") %>% 
  dplyr::select(WBGeneID = wbid, Public_name = SYMBOL, Ontology.name = name, Expected:class)
write.xlsx(fileS6, "../SuppFiles/FileS6_TargetGeneClass_Ontology_Results.xlsx")


