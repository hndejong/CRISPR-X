#!/bin/bash
#LF
#DURGA
#CRISPR-X
WORK_DIR=$1
H1=$2
H2=$3
h1=$4
h2=$5
parent=$6
mismatch=$7

cmd='Rscript /users/lfresard/CRISPR-X/scripts/mutation_reads_multiguide_summary.R ${WORK_DIR} ${H1} ${H2} ${h1} ${h2} ${parent} ${mismatch}'
date >$WORK_DIR/Mutread_hotspot_${H1}-${H2}_${h1}-${h2}_summary.log
echo $cmd >>$WORK_DIR/Mutread_hotspot_${H1}-${H2}_${h1}-${h2}_summary.log
eval $cmd >> $WORK_DIR/Mutread_hotspot_${H1}-${H2}_${h1}-${h2}_summary.log 2>&1
date >>$WORK_DIR/Mutread_hotspot_${H1}-${H2}_summary.log
