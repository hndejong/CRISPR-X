#!bin/bash
#LF
#durga
# January 31st
N1=$1
N2=$2
mismatch=$3
H1=$4
H2=$5
WORK_DIR=$6
BAM_DIR=$7
chr=$8


for i in $(seq {$N1 $N2}); do 
	cmd='bash /users/lfresard/CRISPR-X/scripts/mutations_hotspots.sh ${sample} $mismatch $i $(($i+40)) ${WORK_DIR} ${BAM_DIR} ${chr}'
	echo $cmd >$WORK_DIR/CX${sample}_mutread_hotspot_${i}-${(($i+40))}.log
	eval $cmd >> $WORK_DIR/CX${sample}_mutread_hotspot_${i}-${(($i+40))}.log 2>&1 &

done