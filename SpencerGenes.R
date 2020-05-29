library(dplyr)
library(stringr)
library(UpSetR)

NASdir='/Volumes/onishlab_shared/PROJECTS/15_TF_TEAM/RWC19_Intestine_TF_List/DATASETS/16_Spencer_et_al_2010_FACS_and_pulldown_tilling_array'
if(file.access(NASdir)!=0) { stop("NAS dir does not exist... STOPPING")}

LEpath = paste0(NASdir,'/','LE-intestine_enr_vs_ref.WS200.txt')
L2path = paste0(NASdir,'/','L2-intestine_enr_vs_ref.WS200.txt')

spencerL2genes <-
  read.table(
    L2path,
    quote = "\"",
    comment.char = "",
    header = TRUE
  )
colnames(spencerL2genes) <-
  str_c("spencer_L2_", colnames(spencerL2genes))

spencer_L2_subset <- spencerL2genes %>%
  select(spencer_L2_ID,
         spencer_L2_AveExpr,
         spencer_L2_adj_P_Val,
         spencer_L2_FC)

spencer_L2_subset %>% head

spencerLEgenes <-
  read.table(
    LEpath,
    quote = "\"",
    comment.char = "",
    header = TRUE
  )
colnames(spencerLEgenes) <-
  str_c("spencer_LE_", colnames(spencerLEgenes))

spencer_LE_subset <- spencerLEgenes %>%
  select(spencer_LE_ID,
         spencer_LE_AveExpr,
         spencer_LE_adj_P_Val,
         spencer_LE_FC)

spencer_LE_subset %>% head

attach(spencer_LE_subset)
attach(spencer_L2_subset)

# LEL2_union = union(spencer_L2_ID,spencer_LE_ID)
# LE = LEL2_union %in% spencer_LE_ID
# L2 = LEL2_union %in% spencer_L2_ID
# LEL2 = data.frame(LE=ifelse(LE,1,0),L2=ifelse(L2,1,0))

source('~/work/onish_ChIP_R_Analysis/matchtable.R')
LEL2 = matchtable(spenc_LE=spencer_LE_ID, spenc_L2=spencer_L2_ID)
upset(LEL2)