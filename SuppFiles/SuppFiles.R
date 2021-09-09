library(tidyverse)

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
elt2_peaks_merged <- elt2_peaks %>% left_join(figure3df, by = "WBGeneID") %>% left_join(figure4df, by = "WBGeneID")

write_csv(elt2_peaks_merged, file = "../SuppFiles/FileS1_ELT2_ChIP-seq.csv")

# Give better column names

# FILE S2

fileS2 <- read_csv("./Figure_3_Intestine_Expression_Patterns/02_output/Fig3_Public_Intestine_RNA_Gene_List.csv")
write_csv(fileS2, "../SuppFiles/FileS2_Intestine_Gene_List.csv")

# FILE S3

fileS3 <- read_csv("./Figure_3_Intestine_Expression_Patterns/02_output/FileS3_Intestine-specific_time-resolved_data.csv") %>% 
  dplyr::select(WBGeneID, Sequence.name, Public_name:cluster.description) %>% 
  left_join(figure3df) %>%
  filter(!is.na(intestine_dev_cluster))

write_csv(fileS3, "../SuppFiles/FileS3_Intestine-specific_time-resolved_data.csv")

# FILE S4

fileS4 <- read_csv("./Figure_4_L1_Regulation/02_output/Figure4_WBGeneID_Annotation_Dataframe.csv")
