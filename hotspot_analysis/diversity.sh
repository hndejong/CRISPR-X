MUT_DIR=$1
WORK_DIR=$2
BAM_DIR=$3
sample=$4
mismatch=$5
H1=$6
H2=$7
chr=$8
snp_pos_script=/users/lfresard/CRISPR-X/scripts/snp_position_on_reads_with_allele_identification.pl
summarize_script=/users/lfresard/CRISPR-X/scripts/summarize_diversity_ok.py

cd $WORK_DIR

#1/FILTER BAM FOR READS IN THE HOSTPOT
echo "------------------------------------------"
echo "1/FILTER BAM FOR READS IN THE HOSTPOT"
echo "------------------------------------------"

samtools view -h $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted.bam ${chr}:${H1}-${H2} > $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}.sam
wait

echo "Filter reads spanning hotspot done"
echo "------------------------------------------"


#get list of reads with mutations in hotspot
tail -n+2 $MUT_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}.sam_mutations_per_reads.txt  |awk '$7!=0 {print $1}' > $WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_list.txt


#sort list of read and filetred sam by coord
sort $WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_list.txt >$WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_list_sorted.txt

#samtools sort -n -O sam -o $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_sorted.sam -T $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_sorted_temp $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}.sam 


#intersect them and add the header back
#samtools view -H $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_sorted.sam >$WORK_DIR/header.txt
#cat $WORK_DIR/header.txt > CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot.sam
#grep -v @ $BAM_DIR/CX${n}_n${mismatch}_mapq30_sorted_${H1}-${H2}_sorted.sam >$BAM_DIR/CX${n}_n${mismatch}_mapq30_sorted_${H1}-${H2}_sorted_noheader.sam
#join -j 1 -t $'\t' $BAM_DIR/CX${n}_n${mismatch}_mapq30_sorted_${H1}-${H2}_sorted_noheader.sam $WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_list_sorted.txt >>$WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot.sam




#2/CLEAN BAM FILE
echo "------------------------------------------"
echo "CLEAN BAM FILE"
echo "------------------------------------------"

java -Xmx8g -jar /usr/local/bin/picard.jar CleanSam \
I=$BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}.sam  \
O=$BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_clean.sam 

wait


#3/FILTER SAMREADS WITH LIST OF READS THAT WE WANT
echo "------------------------------------------"
echo "FILTER SAMREADS WITH LIST OF READS THAT WE WANT"
echo "------------------------------------------"

java -Xmx8g -jar /usr/local/bin/picard.jar FilterSamReads \
I=$BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_clean.sam  \
O=$BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_clean_mutreadinhotspot.sam  \
READ_LIST_FILE=$WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_list_sorted.txt \
FILTER=includeReadList
wait

rm $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}.sam


samtools sort $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_clean_mutreadinhotspot.sam $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_clean_mutreadinhotspot_sorted
wait
samtools index $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_clean_mutreadinhotspot_sorted.bam

wait

rm $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_clean.sam
#4/ GET EACH OBSERVED BASE BY READ ON HOTSPOT
echo "------------------------------------------"
echo "GET EACH OBSERVED BASE BY READ ON HOTSPOT"
echo "------------------------------------------"
#position on read.pl
perl $snp_pos_script -v $MUT_DIR/CX${sample}_n${mismatch}_mapq30_sorted_qual30_${H1}-${H2}_filtered.vcf $BAM_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_clean_mutreadinhotspot_sorted.bam >$WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_posfile.txt

#sort by read and position
sort -k2,2 -k1,1 $WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_posfile.txt > $WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_posfile_sorted.txt

wait
rm $WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_posfile.txt
#get only reads that have a mutation somewhere
#cat CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_posfile_sorted.txt |awk '$7!=$8 {OFS="\t"; print $2, $4}'|sort -k1,1 -k2,2 |uniq >CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_posfile_sorted_mutreadsid-direction.txt
wait
#resume for each read what we see
#5/SUMMARIZE RESULTS
echo "------------------------------------------"
echo "SUMMARIZE RESULTS"
echo "------------------------------------------"

python $summarize_script $WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_posfile_sorted.txt

echo "ALL DONE!!!"

wait

rm $WORK_DIR/CX${sample}_n${mismatch}_mapq30_sorted_${H1}-${H2}_mutreadinhotspot_posfile_sorted.txt
