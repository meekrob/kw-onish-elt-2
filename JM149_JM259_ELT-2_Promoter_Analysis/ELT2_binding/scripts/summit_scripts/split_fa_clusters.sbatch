#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:10:00
#SBATCH --qos=normal
#SBATCH --partition=shas
NTHREADS=${SLURM_NTASKS} # passes --ntasks set above
echo "[$0] $SLURM_JOB_NAME $@" # log the command line
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR
date # timestamp

# SLURM_JOB_ID
# SLURM_ARRAY_JOB_ID
# SLURM_ARRAY_TASK_ID is the index variable in array mode

# erinnishgrp@colostate.edu projects setup
PROJ_DIR=/projects/dcking@colostate.edu
source $PROJ_DIR/paths.bashrc

set -x

./split_fa.py allStagesUNION.IDR_0.05.sorted_N_500flank.fa <(cut -f 4 clusters.bed)

cat allStagesUNION.IDR_0.05.sorted_N_500flank_{1..4}.fa > allStagesUNION.IDR_0.05_N_500flank.sorted_not_0.fa
cat allStagesUNION.IDR_0.05.sorted_N_500flank_{0,{2..4}}.fa > allStagesUNION.IDR_0.05.sorted_not_1.fa
cat allStagesUNION.IDR_0.05.sorted_N_500flank_{{0..1},{3..4}}.fa  > allStagesUNION.IDR_0.05.sorted_not_2.fa
cat allStagesUNION.IDR_0.05.sorted_N_500flank_{{0..2},4}.fa > allStagesUNION.IDR_0.05.sorted_not_3.fa
cat allStagesUNION.IDR_0.05.sorted_N_500flank_{0..3}.fa > allStagesUNION.IDR_0.05.sorted_not_4.fa
