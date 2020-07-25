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

# SBATCH ARRAY ARGUMENTS
# ---------------------------
# to do 000,100,200,300,... 
# --array=0-900:100
# To do all 1000 in 0/ with 10 jobs (100 per job):
# sbatch --array=0-900:100 do_compute_features.sbatch 0 100
# To do all 1000 in 1/ with 10 jobs (100 per job):
# sbatch --array=1-1900:100 do_compute_features.sbatch 1 100

# erinnishgrp@colostate.edu projects setup
PROJ_DIR=/projects/dcking@colostate.edu
source $PROJ_DIR/paths.bashrc
module load R/3.5.0

basedir=$1
do_how_many=$2
base_num=$SLURM_ARRAY_TASK_ID
end_num=$(( $base_num + $do_how_many -1 ))

fileargs=$(printf "$basedir/perm_$basedir%03d.bed " $(eval echo {$base_num..$end_num}))
ls -l  $fileargs
outname=$(printf "features-%02d-%03d-%03d.tbl" $basedir $base_num $end_num)

cmd="Rscript compile_geneFeatureMappings.Rscript $outname $fileargs"
echo $cmd
time eval $cmd


