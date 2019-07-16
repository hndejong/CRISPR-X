CRISPR-X BAM files processing
-----------------------------
The scripts reflects analyses performed in the [CRISPR-X](https://www.nature.com/articles/nmeth.4038) paper.


## 1. BAM FILTER ON MAPQ

This first step is filtering bam files on mapping quality, keeping only reads with a mapping quality greater than 30.

Parameters:
* `BAM_DIR`: directory containing BAM files. Naming convention: `CX<sample_number>_n<number_mismatch>.bam` (example: CX10_n5.bam, CX stand for CRISPR-X)
* `FILTERED_BAM_DIR`: output directory containing filtered BAM files
* `REF.fa`: reference genome used in the alignment.
* `SAMPLE_START`: first sample number to start the processing on (works with the current naming of CRISPR-X files, here starting with CX)
* `SAMPLE_STOP`: last sample number on which to perform the analyses
* `mismatches`: number of mismatches allowed in aligned reads

Script: [launch_filterbam.sh](./process_bam/launch_filterbam.sh)
This scripts filters bam files on mapping quality for all samples between `SAMPLE_START` and `SAMPLE_STOP`
```
bash launch_filterbam.sh <BAM_DIR> <FILTERED_BAM_DIR> <REF.fa> <SAMPLE_START> <SAMPLE_STOP> <mismatches>
```
Example command line: `bash launch_filterbam.sh ./BAM ./FILTERED_BAM MYREF.FA 1 20 5`.

This will filter on mapping quality > 30 from sample CX1_n5.bam to CX20_n5.bam and produce file CX1_n5_mapq30_sorted.bam to CX20_n5_mapq30_sorted.bam

This script is using [Filterbam_MAPQ.sh](./process_bam/Filterbam_MAPQ.sh)


## 2. GET ALL VARIATIONS AND COUNTS
This second steps is processing the filtered bams and counting all the alleles observed in the region of interest.
It first generates [mpileup files](http://samtools.sourceforge.net/pileup.shtml). Then it uses another [script](./Get_variations/mpileup_count.py) to generate the allele count observed at each position like in this [example](./count_file_example.md).


Parameters:
* `REF.fa`: reference genome used in the alignment
* `FILTERED_BAM_DIR`: directory containing filtered BAM files
* `OUT_DIR`: output directory 
* `mismatches`:  number of mismatches allowed in aligned reads
* `SAMPLE_START`: first sample number to start the processing on
* `SAMPLE_STOP`: last sample number on which to perform the analyses
* `QUAL`: Base quality used in samtools mpileup (30 used for paper)

Script: [launch_mpileupanalysis.sh](./Get_variations/launch_mpileupanalysis.sh)

```
bash launch_mpileupanalysis.sh <REF.fa> <FILTERED_BAM_DIR> <OUT_DIR> <mismatches> <SAMPLE_START> <SAMPLE_STOP> <QUAL>
```
Example command line: `bash launch_mpileupanalysis.sh MY_REF.fa ./FILTERED_BAM ./COUNTS 5 1 20 30`
This will get allelic counts from the filtered bam files from sample CX1_n5_mapq30_sorted.bam to CX20_n5_mapq30_sorted.bam. Output will look like CX20_n5_mapq30_sorted_qual30.count

This script is using [mpileup_analysis.sh](./Get_variations/mpileup_analysis.sh) which is calling [samtools mpileup](http://www.htslib.org/doc/samtools-1.2.html) and [mpileup_count.py](./Get_variations/mpileup_count.py)


## 3. GET ENRICHMENT PLOT
Parameters:
* `SAMPLE_START`: first sample number to start the processing on
* `SAMPLE_STOP`: last sample number on which to perform the analyses
* `mismatches`:  number of mismatches allowed in aligned reads (5 in [CRISPR-X](https://www.nature.com/articles/nmeth.4038) paper)
* `PARENT_COUNT_FILE`: Reference count file from a sample without CRISPR-X)
* `LOCUS`: chromosome/locus to perform plot on
* `WORK_DIR`: WORKING DIRECTORY

Script: [launch_enrichmentplot.sh](./enrichment/launch_enrichmentplot.sh)
```
bash launch_enrichmentplot.sh <SAMPLE_START> <SAMPLE_STOP> <mismatches> <PARENT_COUNT_FILE> <LOCUS> <WORK_DIR>
```

This script is using [mpileupcounts_plotter_enrichmentonly.R](./enrichment/mpileupcounts_plotter_enrichmentonly.R)

