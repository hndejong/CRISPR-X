#!/bin/bash


REF=$1
BAM_DIR=$2
WORK_DIR=$3
mismatch=$4
N=$5
QUAL=$6
count_SCRIPT=/users/lfresard/CRISPR-X/scripts/mpileup_count.py
DEPTH=1000000

samtools mpileup -d $DEPTH -Q $QUAL -BAf $REF $BAM_DIR/CX${N}_n${mismatch}_mapq30_sorted.bam>$WORK_DIR/CX${N}_n${mismatch}_mapq30_sorted_qual${QUAL}.mpileup
python $count_SCRIPT $WORK_DIR/CX${N}_n${mismatch}_mapq30_sorted_qual${QUAL}.mpileup > $WORK_DIR/CX${N}_n${mismatch}_mapq30_sorted_qual${QUAL}.count
