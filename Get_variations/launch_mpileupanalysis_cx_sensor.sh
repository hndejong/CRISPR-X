#!/bin/bash
BAM_DIR=$1
OUT_DIR=$2
REF=$3



mpileup_script=/users/lfresard/repos/crispr-x/Get_variations/mpileup_analysis_cx_sensor.sh
cd $BAM_DIR

ls *_mapq30_sorted.bam | sed 's/.bam/\t/' | awk '{print $1}' | awk -v BAM_DIR=$BAM_DIR -v OUT_DIR=$OUT_DIR -v REF=$REF 'BEGIN{OFS="\t"}{print $1".bam", OUT_DIR"/"$1".mpileup", OUT_DIR"/"$1".count", REF}' |
	parallel --jobs 10 --col-sep "\t" "${mpileup_script} {1} {2} {3} {4}"