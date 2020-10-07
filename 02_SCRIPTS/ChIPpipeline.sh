#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=0:35:00
#SBATCH --partition=shas
# --open-mode=append assures that the log files will be appended to from different jobs
# These directives will serve as defaults when submitted via sbatch
# but are comments when run via bash
NTHREADS=${SLURM_NTASKS} # passes --ntasks set above
echo "$SLURM_JOB_NAME[$SLURM_JOB_ID] $@" # log the command line
SUBMIT=$0
#JOBSTEPS="BWA BAM SPP IDR UNION"
JOBSTEPS="SPP IDR LOG"
########################################################
##### CONFIGURATION VARIABLES
FLANK=150 # for bedToBw
# the BWA indexes 
BWA_GENOME=/projects/dcking@colostate.edu/support_data/bwa-index/ce11.unmasked.fa
CHROMLENGTHS=/projects/dcking@colostate.edu/support_data/ce11/ce11.chrom.sizes
ALL_STAGES_UNION=allStagesUnion.bed
ALL_STAGES_AGGREGATE=${ALL_STAGES_UNION}.df
# PROJECT ORGANIZATION
INPUT_DIR=01_FASTQ
ALIGN_DIR=03_ALIGN
SPP_DIR=04_SPP
IDR_DIR=05_IDR
SIG_DIR=06_SIGNAL

for dir in $ALIGN_DIR $SPP_DIR $IDR_DIR $SIG_DIR
do
    mkdir -pv $dir
done

# Helper functions
errecho()
{
    1>&2 echo $@
}
run()
{
    echo "running $@"
    time eval $@
}
sb()
{
    sbatch --parsable $@
}
deps() 
{
    vars="$@"
    vars=${vars// /,}
    d=""
    [ -n "$vars" ] && d="-d afterok:$vars"
    echo $d
}
makeFlankFilename() 
{   # anticipate the output filename written by bedToBw.sh 
    # (it is not available as an argument)
    infile=$1
    flank=$2
    root=${infile%.*}
    echo ${root}x${flank}n.bw
}

if [ -z "$SLURM_JOB_ID" ]
######## THIS PART OF THE SCRIPT RUNS IN A NORMAL BASH SESSION    ########
######## but launches JOBSTEPS using sbatch (also in this script) ######## 
######## The job steps are invoked by specifying input and output ######## 
######## filenames only. The specific arguments to the programs   ########
########  called within the jobstep definitions in the sbatch     ########
######## block of this script.                                    ########
then
    # Metadata file
    metadatafile=$1
    shift

    [ $# -gt 0 ] && JOBSTEPS="$@"

    # the "filename" portion of each block must always run 
    # in order to make it possible for the pipeline to be 
    # called with an internal starting point, as when
    # the pipeline fails at some specific step.
    all_ids=""
    idr_filesnames=""
    IDR_JOB_IDS=""
    SCORE_FILEPATHS=""
    SCORE_JOB_IDS="" 
    while read label rep1_fastq rep2_fastq input_fastq
    do
        stage_jids=""
        # skip comment lines
        if [ ${label:0:1} == '#' ] 
        then
            errecho "skipping: $label $rep1_fastq $rep2_fastq $input_fastq" # essentially print 'line' back out
            continue
        fi

        BASE_LOGFILE="$label.jobs.log.$(date +%y%m%d.%H%M)"
        date > $BASE_LOGFILE
        echo -e "\n-----------------------------------------------" >> $BASE_LOGFILE
        echo -e "${label}\t${rep1_fastq}\t${rep2_fastq}\t${input_fastq}" >> $BASE_LOGFILE

        #ALIGN
        # FILENAMES
        rep1_sam=${rep1_fastq/.fastq/.sam}
        rep2_sam=${rep2_fastq/.fastq/.sam}
        input_sam=${input_fastq/.fastq/.sam}
        # JOBS
        if [[ " $JOBSTEPS " =~ " BWA " ]] 
        then
            T="--time=0:05:00"
            align_jid1=$(sb --job-name=bwa-${label}-1 $SUBMIT BWA ${INPUT_DIR}/$rep1_fastq  ${ALIGN_DIR}/$rep1_sam)
            align_jid2=$(sb --job-name=bwa-${label}-2 $SUBMIT BWA ${INPUT_DIR}/$rep2_fastq  ${ALIGN_DIR}/$rep2_sam)
            align_jid3=$(sb --job-name=bwa-${label}-i $SUBMIT BWA ${INPUT_DIR}/$input_fastq ${ALIGN_DIR}/$input_sam)
            stage_jids="$stage_jids $align_jid1 $align_jid2 $align_jid3"
        fi

        #CONVERT FORMAT
        # FILENAMES
        rep1_bam=${rep1_sam/.sam/.bam}
        rep2_bam=${rep2_sam/.sam/.bam}
        input_bam=${input_sam/.sam/.bam}
        # JOBS
        if [[ " $JOBSTEPS " =~ " BAM " ]] 
        then
            D=$(deps $align_jid1)
            bam_jid1=$(sb --job-name=bam-${label}-1 $O $D $SUBMIT BAM $ALIGN_DIR/$rep1_sam $ALIGN_DIR/$rep1_bam)
            D=$(deps $align_jid2)
            bam_jid2=$(sb --job-name=bam-${label}-2 $O $D $SUBMIT BAM $ALIGN_DIR/$rep2_sam $ALIGN_DIR/$rep2_bam)
            D=$(deps $align_jid3)
            bam_jid3=$(sb --job-name=bam-${label}-i $O $D $SUBMIT BAM $ALIGN_DIR/$input_sam $ALIGN_DIR/$input_bam)
            stage_jids="$stage_jids $bam_jid1 $bam_jid2 $bam_jid3"
        fi

        #BW COMPUTE SIGNAL FILES
        # filenames
        # see function: makeFlankFilename()
        rep1_bw=${label}_1.bw
        rep2_bw=${label}_2.bw
        input_bw=${label}_input.bw
        # JOBS
        if [[ " $JOBSTEPS " =~ " BW " ]] 
        then
            ntasks="--ntasks=4"
            tim="--time=0:11:00"
            D=$(deps $bam_jid1)
            bw_jid1=$(sb --job-name=bw-${label}-1 $ntasks $tim $D $SUBMIT BW $ALIGN_DIR/$rep1_bam $SIG_DIR/$rep1_bw)
            D=$(deps $bam_jid2)
            bw_jid2=$(sb --job-name=bw-${label}-2 $ntasks $tim $D $SUBMIT BW $ALIGN_DIR/$rep2_bam $SIG_DIR/$rep2_bw)
            D=$(deps $bam_jid3)
            bw_jid3=$(sb --job-name=bw-${label}-i $ntasks $tim $D $SUBMIT BW $ALIGN_DIR/$input_bam $SIG_DIR/$input_bw)

            stage_jids="$stage_jids $bw_jid1 $bw_jid2 $bw_jid3"
        fi
        
        #BW-SUBTRACT
        # 
        rep1_input_subtracted_bw="${label}_1_minus_input.bw"
        rep2_input_subtracted_bw="${label}_2_minus_input.bw"
        SCORE_FILEPATHS="$SCORE_FILEPATHS $SIG_DIR/$rep1_input_subtracted_bw $SIG_DIR/$rep2_input_subtracted_bw"
        if [[ " $JOBSTEPS " =~ " BW-SUBTRACT " ]] 
        then
            ntasks="--ntasks=4"
            tim="--time=0:03:00"
            D=$(deps $bw_jid1 $bw_jid3)
            bws_jid1=$(sb --job-name=${label}_1-${label}i $ntasks $tim $D $SUBMIT BW-SUBTRACT $SIG_DIR/$rep1_bw $SIG_DIR/$input_bw $SIG_DIR/$rep1_input_subtracted_bw)
            D=$(deps $bw_jid2 $bw_jid3)
            bws_jid2=$(sb --job-name=${label}_2-${label}i $ntasks $tim $D $SUBMIT BW-SUBTRACT $SIG_DIR/$rep2_bw $SIG_DIR/$input_bw $SIG_DIR/$rep2_input_subtracted_bw)

            stage_jids="$stage_jids $bws_jid1 $bws_jid2"
            SCORE_JOB_IDS="$SCORE_JOB_IDS $bws_jid1 $bws_jid2"
        fi

        # PEAK CALLS (SPP)
        # SPP OUTPUT FILENAMES are not specified directly, but have the following format 
        # (like L1_1_VS_L1_input.regionPeak.gz), 
        rep1_regionPeak=${rep1_bam%%.bam}_VS_${input_bam%%.bam}.regionPeak.gz
        rep2_regionPeak=${rep2_bam%%.bam}_VS_${input_bam%%.bam}.regionPeak.gz

        # SPP JOBS
        if [[ " $JOBSTEPS " =~ " SPP " ]] 
        then
            echo "pipeline: $label SPP"
            D=$(deps $bam_jid1 $bam_jid3)
            spp_jid1=$(sb --job-name=spp-${label}-1 $O $D $SUBMIT SPP $label ${ALIGN_DIR}/$rep1_bam ${ALIGN_DIR}/$input_bam $rep1_regionPeak)

            D=$(deps $bam_jid2 $bam_jid3)
            spp_jid2=$(sb --job-name=spp-${label}-2 $O $D $SUBMIT SPP $label ${ALIGN_DIR}/$rep2_bam ${ALIGN_DIR}/$input_bam $rep2_regionPeak)

            stage_jids="$stage_jids $spp_jid1 $spp_jid2"
        fi

        # IDR launch
        # IDR FILENAMES
        idr_out="${label}.narrowPeak"
        idr_filenames="$idr_filenames ${IDR_DIR}/$idr_out"
        # IDR JOBS
        if [[ " $JOBSTEPS " =~ " IDR " ]] 
        then
            echo "pipeline: $label IDR"
            D=$(deps $spp_jid1 $spp_jid2)
            idr_jid=$(sb --ntasks=1 --time=0:02:00 --job-name=idr-${label} $D $SUBMIT IDR ${SPP_DIR}/$rep1_regionPeak ${SPP_DIR}/$rep2_regionPeak ${IDR_DIR}/$idr_out)
            stage_jids="$stage_jids $idr_jid"
        fi

        # Use IDR as the rejoin point from all of the branches.
        IDR_JOB_IDS="$IDR_JOB_IDS $idr_jid"

        #if [[ " $JOBSTEPS " =~ " LOG " ]] 
        if true
        then
            echo "pipeline: $label LOG "
            echo $stage_jids
            # Add a command to merge the temporary log files to the base log file
            echo "# To merge the individual jobs logs into this file:" >> $BASE_LOGFILE
            lognames=""
            for jid in $stage_jids
            do
                lognames="$lognames slurm-$jid.out"
            done
            #D=$(deps $stage_ids)
            #log_jid=$(sb --ntasks=1 --time=0:01:00 --job-name=$label.catlg --output=${label}.catlogs-%j.out $D $SUBMIT LOG $BASE_LOGFILE $lognames)
            if [ -n "$lognames" ]
            then
                echo "cat $lognames >> $BASE_LOGFILE && rm -v $lognames" | tee -a  $BASE_LOGFILE
            fi
        fi

        all_ids="$all_ids $stage_jids"
        
    done < $metadatafile

    # UNION 
    # FILENAMES
    if [[ " $JOBSTEPS " =~ " UNION " ]] 
    then
        D=$(deps $IDR_JOB_IDS)
        union_jid=$(sb $D --ntasks=1 --time=0:05:00 --job-name=union --output=union.%j.out $SUBMIT UNION $ALL_STAGES_UNION  $idr_filenames)
    fi
    # AGGREGATE
    if [[ " $JOBSTEPS " =~ " AGGREGATE " ]]
    then
        D=$(deps $SCORE_JOB_IDS)
        aggregate_jid=$(sb $D  --ntasks=1 --time=0:05:00 --job-name=aggregate $SUBMIT AGGREGATE $ALL_STAGES_UNION $ALL_STAGES_AGGREGATE $SCORE_FILEPATHS)
    fi

    echo "ALL JOBS SUBMITTED:"
    all_ids=$(echo $all_ids) # trim ws
    echo "jid=${all_ids// /,}"

else 
######## THIS PART OF THE SCRIPT RUNS INSIDE SLURM, AND IS CALLED   ########
######## FROM THE BASH SESSION.                                     ########
    export TMPDIR=$SLURM_SCRATCH
    export TMP=$TMPDIR
    jobstep=$1
    shift
    errecho "$jobstep SLURM_JOB_ID=$SLURM_JOB_ID"
    date

    source /projects/dcking@colostate.edu/paths.bashrc
#BWA ###################################
    if [ $jobstep == "BWA" ]
    then
        2>&1 bwa | grep ^Program
        2>&1 bwa | grep ^Version
        infastq=$1
        outsam=$2
        cmd="bwa mem -t $SLURM_NTASKS $BWA_GENOME $infastq > $outsam"
        run $cmd

#BAM ###################################
    elif [ $jobstep == "BAM" ]
    then
        samtools --version
        insam=$1
        outbam=$2
        filter="-F 1536"
        quality="-q 30"
        sort_prefix="$TMPDIR/samsort_$SLURM_JOB_ID"
        cmd="samtools view -@ $SLURM_NTASKS -b $filter $quality -S ${insam} | samtools sort -T $sort_prefix -@ $SLURM_NTASKS -o ${outbam}"
        run $cmd
        cmd="samtools index $outbam"
        run $cmd
        cmd="samtools quickcheck $outbam && rm -v $insam"
        run $cmd

#BW  ###################################
    elif [ $jobstep == "BW" ]
    then
        inbam_path=$1
        inbam_filename=$(basename $inbam_path)
        indir=$(dirname $inbam_path)
        logbw=$2 
        logbw_filename=$(basename $logbw)
        outdir=$(dirname $logbw)
        outbed=${inbam_filename/.bam/.nonlog.bed} 
        # bedToBw.sh - creates a file according to its run parameters
        bamToBw_outfile=$(makeFlankFilename $outbed $FLANK) 
        # example $bamToBw_outfile - LE_nonlogx150n.bw

        # wig file - output of wigmath.LogTransform
        logwig=${bamToBw_outfile/nonlogx${FLANK}n.bw/log.wig}

        # get the chromosome locations from the alignment
        cmd="bedtools bamtobed -i $inbam_path > $outdir/$outbed"
        run $cmd

        # bedToBw.sh- pad the chromosome locations according to $FLANK, scale by read depth.
        # creates $outdir/$bamToBw_outfile
        cmd="bedToBw.sh $outdir/$outbed $FLANK $CHROMLENGTHS -n -bw && rm -v $outdir/$outbed" 
        run $cmd
        
        # perform a log transformation
        # creates a wiggle file (ASCII)
        cmd="meekrob-javaGenomicsToolkit wigmath.LogTransform -p $SLURM_NTASKS -i $outdir/$bamToBw_outfile -o $outdir/$logwig"
        run $cmd

        # convert to bigWig (binary, compressed)
        cmd="wigToBigWig $outdir/$logwig $CHROMLENGTHS $outdir/$logbw_filename && rm -v $outdir/$logwig"
        run $cmd
        
#BW-SUBTRACT ###########################
    elif [ $jobstep == "BW-SUBTRACT" ]
    then
        minuend=$1
        subtrahend=$2
        outfilename=$3
        wigfilename=${outfilename/bw/wig}
        # perform subtraction
        # creates a wiggle file (ASCII)
        cmd="meekrob-javaGenomicsToolkit wigmath.Subtract -z -p $SLURM_NTASKS -m $minuend -s $subtrahend -o $wigfilename"
        run $cmd
        
        # convert to bigWig (binary, compressed)
        cmd="wigToBigWig -clip $wigfilename $CHROMLENGTHS $outfilename && rm -v $wigfilename"
        run $cmd

# AGGREGATE ############################
    elif [ $jobstep == "AGGREGATE" ]
    then
        bedfile=$1
        shift
        outfile=$1
        shift
        scorefiles="$@"
        cmd="meekrob-javaGenomicsToolkit ngs.SplitWigIntervalsToDataFrame -s -l $bedfile -o $outfile $scorefiles"
        run $cmd

#SPP ###################################
    elif [ $jobstep == "SPP" ]
    then
        prefix=$1
        rep=$2
        input=$3
        output=$4
        outdir=${SPP_DIR}
        FRAGLEN=150
        SPP=spp # loadbx has bin/spp as wrapper to run_spp.R
        cmd="$SPP -c=$rep   \
             -i=$input \
             -npeak=300000  \
             -odir=${outdir}  \
             -speak=${FRAGLEN} \
             -p=$SLURM_NTASKS \
             -savr -rf  \
             -savp="${outdir}/$prefix.pdf" \
             -out=${outdir}/$prefix.ccscores"
        run $cmd
        # results do not come out sorted
        sortTempfile="${outdir}/.${output}"
        cmd="zcat $output | sort -k1,1 -k2,2n > $sortTempfile && mv $sortTempfile ${outdir}/$output}"
        run $cmd
        
    
#IDR ###################################
    elif [ $jobstep == "IDR" ]
    then
        rep1=$1
        rep2=$2
        outfile=$3
        outtmp=${outfile/.narrowPeak/.unthresh.narrowPeak}
        IDR_THRESHOLD=0.05
        echo "IDR_THRESHOLD=$IDR_THRESHOLD"
        cmd="idr --samples $rep1 $rep2 --input-file-type narrowPeak \
             --rank signal.value \
             --soft-idr-threshold ${IDR_THRESHOLD} \
             --plot --use-best-multisummit-IDR \
             --random-seed 13 \
             --output-file $outtmp"
        run $cmd
        IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESHOLD} 'BEGIN{print -log(p)/log(10)}')
        echo "IDR_THRESH_TRANSFORMED=$IDR_THRESH_TRANSFORMED"
        idr_filter()
        {
            #awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' $1
            awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $0}' $1
        }
        declare -f idr_filter
        cmd="idr_filter $outtmp | sort -k1 -k2,2n -u > $outfile"
        run $cmd

#BZ ####################################
    elif [ $jobstep == "BW" ]
    then
        errecho "jobstep BW"

#UNION #################################
    elif [ $jobstep == "UNION" ]
    then
        outfile=$1
        shift
        infiles=$@
        cmd="cat $infiles | bed_merge_overlapping.py | sort -k1,1 -k2,2n > $outfile"
        run $cmd
    
#LOG ###################################
    # append the individual slurm logs to the main log file and
    # delete them
    elif [ $jobstep == "LOG" ]
    then
        main=$1
        shift
        outfiles="$@"
        #cmd="cat $outfiles >> $main && rm -v $outfiles"
        cmd="cat $outfiles >> $main"
        run $cmd
#NOT DEFINED
    else
        errecho "jobstep $jobstep is not defined. Must be one of $JOBSTEPS"

    fi # END

fi
