#!/bin/bash
#LF
#Launch Filters on BAM file form a specified directory

filter_script=/users/lfresard/repos/crispr-x/process_bam/Filterbam_MAPQ_cx_sensor.sh
BAM_DIR=$1
OUT_DIR=$2
REF=$3

cd $BAM_DIR

ls *.bam | sed 's/.bam/\t/' | awk '{print $1}' | awk -v BAM_DIR=$BAM_DIR -v OUT_DIR=$OUT_DIR -v REF=$REF 'BEGIN{OFS="\t"}{print $1".bam", OUT_DIR"/"$1"_mapq30.bam", OUT_DIR"/"$1"_mapq30_sorted", REF}' |
	parallel --jobs 10 --col-sep "\t" "${filter_script} {1} {2} {3} {4}"