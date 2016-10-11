#!/bin/bash
mut_n1=$1
mut_n2=$2
mismatch=$3
cont_file=$4
chr=$5
dir=$6

log=CX${mut_n1}_${mut_n2}_${chr}_enrich.log

for i in `seq $mut_n1 $mut_n2` ; do
	cmd="Rscript scripts/mpileupcounts_plotter_enrichmentonly.R CX${i}_n${mismatch}_mapq30_sorted_qual30.count $cont_file $chr $dir"
	echo $cmd >>$dir/$log
	eval $cmd >>$dir/$log 2>&1
done

