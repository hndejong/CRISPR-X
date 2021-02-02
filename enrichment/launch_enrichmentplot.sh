#!/bin/bash

DIR=$1
mut_n1=$2
mut_n2=$3
mismatch=$4
cont_file=$5
chr=$6
dir=$7
count_dir=$8

log=CX${mut_n1}_${mut_n2}_${chr}_enrich.log

for i in `seq $mut_n1 $mut_n2` ; do
	cmd="Rscript $DIR/enrichment/mpileupcounts_plotter_enrichmentonly.R $count_dir/CX${i}_n${mismatch}_mapq30_sorted_qual30.count $cont_file $chr $dir"
	echo $cmd >>$dir/$log
	eval $cmd >>$dir/$log 2>&1
done
