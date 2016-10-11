#!bin/bash
#LF
N1=$1
N2=$2
mismatch=$3
H1=$4
H2=$5
WORK_DIR=$6
BAM_DIR=$7
MUT_DIR=$8
chr=$9

for i in $(seq $N1 $N2) ; do 
	cmd='bash /users/lfresard/CRISPR-X/scripts/diversity.sh ${MUT_DIR} ${WORK_DIR} ${BAM_DIR} ${i} ${mismatch} ${H1} ${H2} ${chr}'
	date >$WORK_DIR/CX${i}_diversity_${H1}-${H2}.log
	echo $cmd >>$WORK_DIR/CX${i}_diversity_${H1}-${H2}.log
	eval $cmd >> $WORK_DIR/CX${i}_diversity_${H1}-${H2}.log 2>&1 &
	date >> $WORK_DIR/CX${i}_diversity_${H1}-${H2}.log
done