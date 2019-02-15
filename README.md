#CRISPR-X BAM files processing.



## 1. BAM FILTER ON MAPQ


```bash launch_filterbam.sh <BAM_DIR> <FILTERED_BAM_DIR> <REF.fa> <SAMPLE_START> <SAMPLE_STOP> <mismatches>```

--> IS CALLING Filterbam_MAPQ.sh
---------------------

2/GET ALL VARIATIONS AND COUNTS
bash launch_mpileupanalysis.sh <REF.fa> <FILTERED_BAM_DIR> <OUT_DIR> <mismatches> <SAMPLE_START> <SAMPLE_STOP> <QUAL>

--> IS CALLING mpileup_analysis.sh THAT IS CALLING samtools mpileup AND mpileup_count.py

---------------------

3/GET ENRICHMENT PLOT
bash launch_enrichmentplot.sh <SAMPLE_START> <SAMPLE_STOP> <mismatches> <PARENT_COUNT_FILE> <LOCUS> <WORK_DIR>
--> IS CALLING mpileupcounts_plotter_enrichmentonly.R

---------------------

SCRIPTS ARE IN DURGA: /users/lfresard/CRISPR-X/scripts

EXAMPLES
1/bash scripts/launch_filterbam.sh /mnt/lab_data/bassik/gaelenh2/CRISPR-X/Results/gaelenh2_CX179-198 /users/lfresard/CRISPR-X/data/Filtered_bam/CX179-198 /users/lfresard/CRISPR-X/data/Reference/GFPCherry684.fa 179 198 5
2/bash scripts/launch_mpileupanalysis.sh /users/lfresard/CRISPR-X/data/Reference/GFPCherry684.fa /users/lfresard/CRISPR-X/data/Filtered_bam/CX179-198 /users/lfresard/CRISPR-X/data/mpileup/CX179-198 5 179 198 30
3/bash scripts/launch_enrichmentplot.sh 179 187 5 CX188_n5_mapq30_sorted_qual30.count GFP684 /users/lfresard/CRISPR-X/data/mpileup/CX179-198
