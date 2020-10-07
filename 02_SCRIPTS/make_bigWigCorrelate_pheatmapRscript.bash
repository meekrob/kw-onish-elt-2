#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --partition=shas
#SBATCH --output=pheatmap-%j.log
#SBATCH --job-name=pheatmap

fileargs=( $@ )
if (( ${#fileargs[@]} < 1))
then
    echo "# Usage: $0 file1.bw file2.bw ..."
    echo "# This script generates an R script with data inserted from all"
    echo "# pairwise comparisons of the given bigwig files using bigWigCorrelate."
    echo "# Required: bigwigCorrelate (UCSC User Apps). R libraries pheatmap, stringr"

    exit
fi

source /projects/dcking@colostate.edu/paths.bashrc

echo "###### BEGIN R SCRIPT #######"
echo "# Command \"$0 $@\""

echo "require(stringr)"
echo "require(pheatmap)"


# go by row, the first element is omitted 
for ((i=1; $i < (( ${#fileargs[@]} )); i++))
do
    if (( i == 1 ))
    then 
        echo -ne "k=c($output"
    fi
    lhf=${fileargs[$i]}

    let j=$i

    # comparison end one sooner each time
    for ((j=0; $j < $i; j++))
    do
        rhf=${fileargs[$j]}
        cmd="bigWigCorrelate -ignoreMissing $lhf $rhf" 
        output=$(eval $cmd)
        if (( i > 1 ))
        then
            echo -ne ", "
        fi

        echo -ne "$output"
    done
done
    echo ")"

echo "# load values as an upper triangular matrix"
n=${#fileargs[@]}
echo "kmx = matrix(1, $n,$n)"
echo "kmx[upper.tri(kmx, diag = F)] <- k"
echo "kmx = kmx * t(kmx) # use the identity in the lower triangle to reflect the upper triangle"
rcnames=$(echo ${fileargs[@]})
echo "rcnames=str_split(\"$rcnames\", ' ')[[1]]"
echo "rownames(kmx)<-rcnames"
echo "colnames(kmx)<-rcnames"
echo "pheatmap(kmx)"

echo "###### END R SCRIPT #######"
