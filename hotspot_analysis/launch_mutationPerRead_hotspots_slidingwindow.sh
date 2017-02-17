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
window_size=$8
mode=$9

for i in $(seq $H1 $H2); do 
	j=$(($i+$window_size))
	cmd='bash /users/lfresard/repos/crispr-x/hotspot_analysis/mutations_hotspots_slidingwindow.sh ${sample} $mismatch $i $j ${WORK_DIR} ${BAM_DIR} ${chr} ${mode}'
	echo $cmd >$WORK_DIR/CX${sample}_mutread_hotspot_${i}-${j}.log
	eval $cmd >> $WORK_DIR/CX${sample}_mutread_hotspot_${i}-${j}.log 2>&1

done