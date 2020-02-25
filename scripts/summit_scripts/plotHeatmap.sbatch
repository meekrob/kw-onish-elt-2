#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --partition=shas
NTHREADS=${SLURM_NTASKS} # passes --ntasks set above
echo "[$0] $SLURM_JOB_NAME $@" # log the command line
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR
date # timestamp

# erinnishgrp@colostate.edu projects setup
PROJ_DIR=/projects/dcking@colostate.edu
source $PROJ_DIR/paths.bashrc

matrixfile=$1

cmd="plotHeatmap -m $matrixfile -out ${matrix_file}_${SLURM_JOB_ID}.png"
echo $cmd
time eval $cmd
