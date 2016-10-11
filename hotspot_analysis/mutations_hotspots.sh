#!/bin/bash
#LF
#Mai 2016
n=$1
mismatch=$2
H1=$3
H2=$4
WORK_DIR=$5
BAM_DIR=$6
chr=$7

date

#1/FILTER BAM FOR READS IN THE HOSTPOT
echo "------------------------------------------"
echo "1/FILTER BAM FOR READS IN THE HOSTPOT"
echo "------------------------------------------"

samtools view $BAM_DIR/CX${n}_n${mismatch}_mapq30_sorted.bam ${chr}:${H1}-${H2} > $BAM_DIR/CX${n}_n${mismatch}_mapq30_sorted_${H1}-${H2}.sam
wait

echo "Filter reads spanning hotspot done"


#2/PROCESS THIS FILTERED BAM INTO PYTHON SCRIPT TO COUNT NUMBER OF MUTATIONS
echo "------------------------------------------"
echo "PROCESS THIS FILTERED SAM INTO PYTHON SCRIPT TO COUNT NUMBER OF MUTATIONS"
echo "------------------------------------------"

cd $WORK_DIR
python /users/lfresard/CRISPR-X/scripts/getnumber_mutations_read.py $BAM_DIR/CX${n}_n${mismatch}_mapq30_sorted_${H1}-${H2}.sam $H1 $H2

echo "mutation counts on reads spanning hotspot done"

date