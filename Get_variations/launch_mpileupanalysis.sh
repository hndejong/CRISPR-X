#!/bin/bash

REF=$1
BAM_DIR=$2
OUT_DIR=$3
mismatch=$4
N1=$5
N2=$6
QUAL=$7


SCRIPT_DIR=/users/lfresard/CRISPR-X/scripts

for i in `seq $N1 $N2`; do
	bash $SCRIPT_DIR/mpileup_analysis.sh $REF $BAM_DIR $OUT_DIR $mismatch $i $QUAL &
done