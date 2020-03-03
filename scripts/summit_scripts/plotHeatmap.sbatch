#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:02:00
#SBATCH --qos=normal
#SBATCH --partition=shas
NTHREADS=${SLURM_NTASKS} # passes --ntasks set above
echo "[$0] $SLURM_JOB_NAME $@" # log the command line
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR
date # timestamp

ml purge
ml singularity
plotHeatmap() {
    singularity exec /projects/dcking@colostate.edu/containers/python3.img plotHeatmap $@
}

matrixfile=$1

#### plotHeatmap
# clear these in the plot, the defaults are TSS and TES
startlabel="--startLabel=''" 
endlabel="--endLabel=''" 

pdffile=${matrixfile/.mx/${SLURM_JOB_ID}.pdf}

categories=LE_specific,Post_embryonic,Dev_increasing,L3_specific,Not_changing
categories_label=$(echo $categories | sed 's/,/ /g')
regionsLabel="--regionsLabel $categories_label"
samplesLabel="--samplesLabel late_embryo larva_stage_1 larva_stage_3"
#cmd="plotHeatmap -m $matrixfile -out ${matrix_file}_${SLURM_JOB_ID}.png"
cmd="plotHeatmap -m $matrixfile $startlabel $endlabel $samplesLabel $regionsLabel -out $pdffile --colorMap RdYlBu_r --heatmapHeight 11 --heatmapWidth 2"
echo $cmd
time eval $cmd
