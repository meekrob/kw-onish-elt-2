for f in $(ls ../03_ALIGN/*.bam); do ln -sv  $(readlink $f) $(basename $f); done
#‘L1_1.bam’ -> ‘/pl/active/onishimura_lab/PROJECTS/04_elt2_ChIPseq_Erin_David/00_modERN_DATA/02_L1/L1_1_ENCFF961UHP.bam’
#‘L1_2.bam’ -> ‘/pl/active/onishimura_lab/PROJECTS/04_elt2_ChIPseq_Erin_David/00_modERN_DATA/02_L1/L1_2_ENCFF531MIE.bam’
#‘L1_input.bam’ -> ‘/pl/active/onishimura_lab/PROJECTS/04_elt2_ChIPseq_Erin_David/00_modERN_DATA/02_L1/L1_input_ENCFF291EVK.bam’
#‘L3_1.bam’ -> ‘/pl/active/onishimura_lab/PROJECTS/04_elt2_ChIPseq_Erin_David/00_modERN_DATA/03_L3/L3_1_ENCFF032LTQ.bam’
#‘L3_2.bam’ -> ‘/pl/active/onishimura_lab/PROJECTS/04_elt2_ChIPseq_Erin_David/00_modERN_DATA/03_L3/L3_2_ENCFF286GQA.bam’
#‘L3_input.bam’ -> ‘/pl/active/onishimura_lab/PROJECTS/04_elt2_ChIPseq_Erin_David/00_modERN_DATA/03_L3/L3_input_ENCFF813CFC.bam’
#‘LE_1.bam’ -> ‘/pl/active/onishimura_lab/PROJECTS/04_elt2_ChIPseq_Erin_David/00_modERN_DATA/01_LE/LE_1_ENCFF716XRG.bam’
#‘LE_2.bam’ -> ‘/pl/active/onishimura_lab/PROJECTS/04_elt2_ChIPseq_Erin_David/00_modERN_DATA/01_LE/LE_2_ENCFF280JAO.bam’
#‘LE_input.bam’ -> ‘/pl/active/onishimura_lab/PROJECTS/04_elt2_ChIPseq_Erin_David/00_modERN_DATA/01_LE/LE_input_ENCFF113NFX.bam’
