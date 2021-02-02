#!/bin/bash

DIR=$1
REF=$2
BAM_DIR=$3
OUT_DIR=$4
mismatch=$5
N1=$6
N2=$7
QUAL=$8

for i in `seq $N1 $N2`; 
do
	bash $DIR/Get_variations/mpileup_analysis.sh $DIR $REF $BAM_DIR $OUT_DIR $mismatch $i $QUAL &
done
