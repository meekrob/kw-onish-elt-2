#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=0:05:00
#SBATCH --qos=normal
#SBATCH --partition=shas
NTHREADS=${SLURM_NTASKS} # passes --ntasks set above
echo "[$0] $SLURM_JOB_NAME $@" # log the command line
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR
date # timestamp

# erinnishgrp@colostate.edu projects setup
PROJ_DIR=/projects/dcking@colostate.edu
#source $PROJ_DIR/paths.bashrc
ml purge
ml singularity
plotHeatmap() {
    singularity exec $PROJ_DIR/containers/python3.img plotHeatmap $@
}

computeMatrix() {
    singularity exec $PROJ_DIR/containers/python3.img computeMatrix $@
}

plotProfile() {
    singularity exec $PROJ_DIR/containers/python3.img plotProfile $@
}

# these will be in the filenames somewhere (see below)
categories=LE_specific,Post_embryonic,Dev_increasing,L3_specific,Not_changing
categories_label=$(echo $categories | sed 's/,/ /g')

#cd threshold_range

# grab all the thresholds in the list of files
#for T in 25 33 50 66 75;
#for T in $(ls allStagesUNION_max_sd*_Not_changing.bedlike | cut -d '_' -f4)
for T in 0.010;
do

    regionfiles=$(eval "echo allStagesUNION_max_sd_${T}_{$categories}.bedlike")
    #regionfiles=peakgenes.WS271.k{0..4}.bed
    matrixfile=${T}_${SLURM_JOB_ID}.mx

    #### computeMatrix
    #scorefiles="../log_LE_1_minus_log_LE_input_z.bw ../log_LE_2_minus_log_LE_input_z.bw ../log_L1_1_minus_log_L1_input_z.bw ../log_L1_2_minus_log_L1_input_z.bw ../log_L3_1_minus_log_L3_input_z.bw ../log_L3_2_minus_log_L3_input_z.bw"
    scorefiles="../LE_signal.bw ../L1_signal.bw ../L3_signal.bw"
    #cmd="computeMatrix scale-regions -p $SLURM_NTASKS -S $scorefiles -R $regionfiles --samplesLabel LE_1 LE_2 L1_1 L1_2 L3_1 L3_2 -o $matrixfile"
    cmd="computeMatrix scale-regions -p $SLURM_NTASKS -a 1000 -b 1000 -S $scorefiles -R $regionfiles --samplesLabel late_embryo larva_stage_1 larva_stage_3 -o $matrixfile"
    echo $cmd
    [ -n "$SLURM_NTASKS" ] && time eval $cmd

    #### plotHeatmap
    startlabel="--startLabel=''" # defaults to TSS
    endlabel="--endLabel=''" # defaults to TES

    pngfile=${matrixfile/.mx/.pdf}
    cmd="plotHeatmap -m $matrixfile $startlabel $endlabel --samplesLabel LE_1 LE_2 L1_1 L1_2 L3_1 L3_2 --regionsLabel $categories_label -out $pngfile"
    cmd="plotHeatmap -m $matrixfile $startlabel $endlabel --samplesLabel late_embryo larva_stage_1 larva_stage_3 --regionsLabel $categories_label -out $pngfile --colorMap RdYlBu_r --heatmapHeight 11 --heatmapWidth 2"
    echo $cmd
    [ -n "$SLURM_NTASKS" ] && time eval $cmd

    #cmd="plotProfile -m $matrixfile -out ExampleProfile1.png "


done
