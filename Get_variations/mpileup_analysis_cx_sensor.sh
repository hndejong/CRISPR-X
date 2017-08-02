#!/bin/bash
# LF
# durga
# august 2017
# 1/perform mpileup analyis on filtered bam
# 2/get allele counts per position from the mpileup file

# Assign variables
BAM_FILE=$1
MPILEUP_FILE=$2
COUNT_FILE=$3
REF=$4

QUAL=30
count_SCRIPT=/users/lfresard/repos/crispr-x/Get_variations/mpileup_count.py
DEPTH=1000000

date

echo "---------"
echo $BAM_FILE

echo "converting bam to mpileup file"
samtools mpileup -d $DEPTH -Q $QUAL -BAf $REF $BAM_FILE >$MPILEUP_FILE

echo "converting mpiluep to count file"
python $count_SCRIPT $MPILEUP_FILE > $COUNT_FILE

echo "---------"

