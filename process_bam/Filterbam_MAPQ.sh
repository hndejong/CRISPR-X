#!/bin/bash
#LF
#Filter for MAPQ>30. 
#Input example: CX1_n5.bam
# Output name example: CX1_n5_mapq30_sorted.bam. This will be the input for samtools mpileup in the next step.

date

WORK_DIR=$1
OUT_DIR=$2
REF=$3

N=$4
mismatch=$5
bam_file=CX${N}_n${mismatch}.bam
mapq_suf=$(basename $bam_file .bam)_mapq30.bam
mapq_sorted=$(basename $mapq_suf .bam)_sorted

cd $OUT_DIR
#1/ MAPQ30
samtools view -h $WORK_DIR/$bam_file | awk '($5>=30)||($1=="@SQ")'| samtools view -Sbh - > $OUT_DIR/$mapq_suf #filter for mapping quality 
samtools sort -m 10G $mapq_suf $mapq_sorted # sort filtered bam
samtools index $mapq_sorted.bam # index filtered bam

cd -
