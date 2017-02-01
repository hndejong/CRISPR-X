#!/bin/bash
#LF
#DURGA
#CRISPR-X
WORK_DIR=$1
parent=$2
sample=$3

cmd='Rscript /users/lfresard/repos/crispr-x/hotspot_analysis/mutation_reads_summary_slidingwindow.R ${WORK_DIR} ${parent} ${sample}'
date >$WORK_DIR/Mutread_hotspot_CX${sample}_slidingwindow_summary.log
echo $cmd >>$WORK_DIR/Mutread_hotspot_CX${sample}_slidingwindow_summary.log
eval $cmd >> $WORK_DIR/Mutread_hotspot_CX${sample}_slidignwindow_summary.log 2>&1
date >>$WORK_DIR/Mutread_hotspot_CX${sample}_slidingwindow_summary.log
