#!/bin/bash
#LF
#1/Filter for MAPQ>30
#2/Filter for PCR duplicates

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
samtools view -h $WORK_DIR/$bam_file | awk '($5>=30)||($1=="@SQ")'| samtools view -Sbh - > $OUT_DIR/$mapq_suf
samtools sort -m 10G $mapq_suf $mapq_sorted
samtools index $mapq_sorted.bam

cd -
