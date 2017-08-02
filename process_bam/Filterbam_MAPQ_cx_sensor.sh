#!/bin/bash
#LF
#August 2017
#1/Filter BAM files for MAPQ>30

date

BAM_FILE=$1
FILTERED_BAM_FILE=$2
SORTED_FILTERED_BAM_FILE=$3
REF=$4


echo "---------"
echo $BAM_FILE
echo "filtering input bam file for mapping quality"
#1/ MAPQ30
samtools view -h $BAM_FILE | awk '($5>=30)||($1=="@SQ")'| samtools view -Sbh - > $FILTERED_BAM_FILE

echo "sorting output by coordinates"
samtools sort -m 10G $FILTERED_BAM_FILE $SORTED_FILTERED_BAM_FILE
samtools index $SORTED_FILTERED_BAM_FILE.bam

echo "---------"

