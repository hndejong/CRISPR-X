#!/bin/bash
#LF
#Launch Filters on Bam form a specified directory

SCRIPT_DIR=/users/lfresard/CRISPR-X/scripts
BAM_DIR=$1
OUT_DIR=$2
REF=$3
n1=$4
n2=$5
mismatch=$6
for i in $(seq $n1 $n2) ; do
	cmd="bash $SCRIPT_DIR/Filterbam_MAPQ.sh $BAM_DIR $OUT_DIR $REF ${i} ${mismatch}"
	echo $cmd
	eval $cmd &
	
done

