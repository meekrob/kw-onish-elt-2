#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=0:25:00
#SBATCH --partition=shas
echo "Using ${SLURM_NTASKS:=1} cores."  # set to 1 if unset (happens in jupyter, sinteractive
echo "$SLURM_JOB_NAME[$SLURM_JOB_ID] $@" # log the command line

bamfile=$1
subsamps="30 35 40 45 50 55 60 65 70 75 80 85 90 95"
seed=1
for samp in $subsamps
do
    echo "Subsample the input alignment ($bamfile) to $samp%:"

    outfile=${bamfile/.bam/.${seed}_${samp}.bam}
    cmd="samtools view -@ $SLURM_NTASKS -s ${seed}.${samp} $bamfile -b -o $outfile"
    echo $cmd
    [ -n "$SLURM_JOB_ID" ] && time eval $cmd # only execute when sbatch'd

    echo
    echo "Count the number subsampled:"
    cmd="nreads=\$(samtools view -@ $SLURM_NTASKS -c $outfile)"
    echo $cmd
    [ -n "$SLURM_JOB_ID" ] && time eval $cmd # only execute when sbatch'd

    echo "BAM $bamfile $outfile $samp $nreads"
    echo

    seed=$(( seed + 1 ))
done
