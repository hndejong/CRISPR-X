#!bin/bash
#LF
N1=$1
N2=$2
mismatch=$3
H1=$4
H2=$5
WORK_DIR=$6
BAM_DIR=$7
chr=$8

for i in $(seq $N1 $N2) ; do 
	cmd='bash /users/lfresard/CRISPR-X/scripts/mutations_hotspots.sh ${i} ${mismatch} ${H1} ${H2} ${WORK_DIR} ${BAM_DIR} ${chr}'
	echo $cmd >$WORK_DIR/CX${i}_mutread_hotspot_${H1}-${H2}.log
	eval $cmd >> $WORK_DIR/CX${i}_mutread_hotspot_${H1}-${H2}.log 2>&1 &

done