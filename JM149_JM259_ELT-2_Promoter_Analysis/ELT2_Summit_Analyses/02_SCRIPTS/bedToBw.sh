#!/usr/bin/env bash

#
# Copyright (c) 2018, David C. King 
# Copyright (c) 2015, Erin Nishimura
#
# All rights reserved.
# Adapted from Kohta's 2012 script: bedToBw.sh 

VERSION=2

printf "bedToBw:   START\n"
DEP_NAMES=(bedtools bedGraphToBigWig)
declare -A DEP_CHECK

check_dependencies() {
    for dep in ${DEP_NAMES[@]}
    do
        if hash $dep
        then
            DEP_CHECK[$dep]="OK"
        else
            DEP_CHECK[$dep]="not found"
        fi

    done
}

if [ $# -lt 5 ]
then
	echo -e "\t`basename $0` version $VERSION"
	echo -e "\tby Erin Nishimura; all rights reserved"
  	echo -e "\n\tUsage: `basename $0` <file.bed> <extension> <Genome> -c/n -r/bg/bw"
	echo -e "\n\t<hits.bed>\tBowtie output. Quality is not considered"
	echo -e "\t<extension>\tDesired final fragment size in bp."
	echo -e "\t<genome>\tSorted chromosome length file."
	echo -e "\t-c/n\t\tBase count (-c) or coverage-normalized base count (-n)."
	echo -e "\t-r/bg/bw\t\tOutput is read length only or bedgraph."
	echo -e "\n\tstdout not available\n"
	echo -e "\n\tRequires bedtools to be preinstalled."
    hash bedtools || echo ".........bedtools is NOT installed."
	echo -e "\tRequires UCSC tools to be preinstalled.\n"
    hash bedGraphToBigWig || echo "UCSC tools is NOT installed."

    check_dependencies
    echo ">>>>>>>>>> ${DEP_CHECK[@]} <<<<<<<<<<"
	exit 
fi

# enforce dependencies
hash bedtools || { echo "bedtools is NOT installed."; exit 1; }
hash bedGraphToBigWig || { echo "bedGraphToBigWig is NOT installed."; exit 1; }
DEBUG_SCRIPTNAME=$(basename $0)
# ANSI escape codes to be used with echo -e
RED='\033[0;31m'
BLACK='\033[0;30m'

# USAGE:
# $0 function name
# $1 hits.bed file in which each entry is a mapped read
# $2 Extension (final length)
# $3 Genome
# $4 Normalize -c/n
# $5 Output format: -r Print only read length; -bg bedGraph; -bw bigwig




# 1) Determine output file names
	F=$1
	seed="${F%%.bed}"	
	if [ $4 = "-c" ]; then	
        echo -e "$RED$DEBUG_SCRIPTNAME: using -c$BLACK"
		bgout="$seed"x"$2"c.bg
		bwout="$seed"x"$2"c.bw
	elif [ $4 = "-n" ]; then	
        echo -e "$RED$DEBUG_SCRIPTNAME: using -n$BLACK"
		bgout="$seed"x"$2"n.bg
		bwout="$seed"x"$2"n.bw
	else 
		echo -e "\tNormalization option (-n/c) required. Exit."
		exit
	fi
        
        #echo $seed
        #printf "\n"
	
        
# 2) Get read length. Get this from the wc of the input bed file ($1)
	rlength=$(head -n 1 $1 | awk '{print($3-$2)}') 
    if [ ! $? ]
    then
        echo -e "$RED$DEBUG_SCRIPT: there was a problem, failing$BLACK"
        exit 1
    
    fi
	
	if [ $5 = "-r" ]; then
		echo -e "bedToBw:   Here is the read length:" $rlength
		exit
	fi
        
        printf "bedToBw:   Read length is $rlength \n"
        
        
# 3) Get the genome length. Get this from the sum of the chromosome length file ($3)
        
        printf "bedToBw:   Getting genome size...\n"
	
        genome=$(cat $3 | awk '{total+=$2}END{print total}')

        echo -e "bedToBw:   Genome size is: $genome"

# 4) Calculate scale factor. If $4 is "-c", scale=1. If $4 is "-n", scale factor is GenomeSize / (final#BpExtendedReads * #ReadsInInputBedFile)
	if [ $4 = "-c" ]; then
		scale=1
	elif [ $4 = "-n" ]; then
                echo "bedToBw:   Getting the number of mapped reads..."
                mappedreads=$(wc -l $1 | awk '{ print $1}')
                scale=$(wc -l $1 | awk '{printf "%.2f", '$genome'/('$2'*$1)}')
                echo -e "bedToBw:   The number of mapped reads is: ${mappedreads}."
	else
		echo -e "\tNormalization option (-n/c) required. Exit."
		exit
	fi
        
        echo -e "bedToBw:   Normalization factor is: $scale"

# 5) Sort using "linux sort"
# 6) Extend reads by final length minus the read length ($2 - $rlength)
# 7) Convert to a bedgraph using "genomeCoverageBed" (bedtools)
        final=$2
        n="$((final - rlength))"
        #printf "bedToBw:   Bed file extend by $n bases.\n"
        
        printf "bedToBw:   Sorting\n"
        tmp_sort_bed=$(mktemp -t bed_sort.XXXXX)
        echo "sort tmp_file is $tmp_sort_bed"

        cmd1="sort -k1,1 -k2,2n $1 > $tmp_sort_bed"
        echo -e "\t\t$cmd1"
        if ! eval time $cmd1
        then 
            echo "FAILED: $cmd1" >&2
            exit 2
        fi
            
        
        printf "bedToBw:   Extending reads by ${n} bases\n"
        #tmp_slop_bed="slop_temp_$1"
        tmp_slop_bed=$(mktemp -t bed_splot.XXXXX)
        echo "slop tmp_file is $tmp_slop_bed"
        cmd2="bedtools slop -s -i $tmp_sort_bed -g ${3} -l 0 -r ${n} > $tmp_slop_bed"
        echo -e "\t\t$cmd2"
        if ! eval time $cmd2
        then
            echo "FAILED: $cmd2" >&2
            exit 2
        fi
        
        printf "bedToBw:   Converting to bedgraph\n"
        cmd3="bedtools genomecov -bg -i $tmp_slop_bed -g ${3} -scale ${scale} > ${bgout}"
        echo -e "\t\t$cmd3"
        if ! eval time $cmd3
        then 
            echo "FAILED: $cmd3" >&2
            exit 2
        fi
        

# 8) Convert bedgraph to a bigwig with bedGraphToBigWig:
        printf "bedToBw:   Converting to bigwig\n"
        
        cmd4="bedGraphToBigWig ${bgout} ${3} ${bwout}"
        echo -e "\t\t$cmd4"
        if ! eval time $cmd4
        then
            echo "FAILED: $cmd3" >&2
            exit 2
        fi

# 9) cleanup ${bgout}
        printf "bedToBw:   Removing temp files\n"
        rm $tmp_sort_bed
        rm $tmp_slop_bed
        rm ${bgout}
        
        printf "bedToBw:   END\n"
        


