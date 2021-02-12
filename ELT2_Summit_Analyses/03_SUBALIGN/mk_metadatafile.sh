echo -n "#$0 "
date
samps="30 35 40 45 50 55 60 65 70 75 80 85 90 95"
for stage in LE L1 L3
do
    for samp in $samps
    do
        rep1=$(grep -h ^BAM *.out | grep ${stage}_1 | grep -w $samp | awk '{print $3}')
        rep2=$(grep -h ^BAM *.out | grep ${stage}_2 | grep -w $samp | awk '{print $3}')
        input=${stage}_input.bam
        echo "${stage}_$samp	${rep1/.bam/.fastq}	${rep2/.bam/.fastq}	${input/.bam/.fastq}"
    done
done
