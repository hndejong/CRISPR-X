#!bin/bash
#LF
#durga
# January 31st
sample=$1
mismatch=$2
H1=$3
H2=$4
WORK_DIR=$5
BAM_DIR=$6
chr=$7


for i in $(seq {$H1 $H2}); do 
	cmd='bash /users/lfresard/CRISPR-X/scripts/mutations_hotspots.sh ${sample} $mismatch $i $(($i+40)) ${WORK_DIR} ${BAM_DIR} ${chr}'
	echo $cmd >$WORK_DIR/CX${sample}_mutread_hotspot_${i}-${(($i+40))}.log
	eval $cmd >> $WORK_DIR/CX${sample}_mutread_hotspot_${i}-${(($i+40))}.log 2>&1 &

done