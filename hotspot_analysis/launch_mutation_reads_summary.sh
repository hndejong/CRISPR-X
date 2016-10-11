#!/bin/bash
#LF
#DURGA
#CRISPR-X
WORK_DIR=$1
H1=$2
H2=$3
parent=$4
mismatch=$5

cmd='Rscript /users/lfresard/CRISPR-X/scripts/mutation_reads_summary.R ${WORK_DIR} ${H1} ${H2} ${parent} ${mismatch}'
date >$WORK_DIR/Mutread_hotspot_${H1}-${H2}_summary.log
echo $cmd >>$WORK_DIR/Mutread_hotspot_${H1}-${H2}_summary.log
eval $cmd >> $WORK_DIR/Mutread_hotspot_${H1}-${H2}_summary.log 2>&1
date >>$WORK_DIR/Mutread_hotspot_${H1}-${H2}_summary.log
