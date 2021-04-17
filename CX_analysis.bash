#!/bin/bash

# Wrapper script written by Hannah De Jong for Ashley Lab Use on Sherlock server, October 2020
# Requires slurm
# Requires prior installation of R packages- see example commented commands
# May require modifications for use on other servers

# Example call:
# CX_analysis.bash data_path ref_file chrom_name code_path number_of_samples ctrl_sample_number

# Requires specific data input naming and path configuration; see "CX_analysis_conventions.txt"

# Store filepaths
DATA=$1 # Path to data
REF_NAME=$2 # Name of reference sequence fasta file
CHROM_NAME=$3 # Name of chromosome for alignment within reference sequence
CODE=$4 # Path to CRISPR-X code
FIRST_SAMPLE=$5 # First sample number
LAST_SAMPLE=$6 # Last sample number 
CTRL_SAMPLE=$7 # File number of control sample

# Create output directories

mkdir -p $DATA/counts
mkdir -p $DATA/enrichment

# Load modules

ml biology
ml samtools
ml python
ml R

# Do alignment
# May want to specify slurm queues with "-p" argument
sbatch --wait $CODE/align_all_paired.bash $DATA $REF_NAME $FIRST_SAMPLE $LAST_SAMPLE

# Make count and pileup files for each sample

for i in $(seq $FIRST_SAMPLE $LAST_SAMPLE)
do
	sbatch --wait $CODE/Get_variations/mpileup_analysis.sh $CODE $DATA/seqs/ref/$REF_NAME $DATA/seqs/filtered_bam $DATA/counts 5 ${i} 30
done

#Install R packages before enrichment script called

#R
#install.packages('RColorBrewer', repos='http://cran.us.r-project.org')
#install.packages('ggplot2', repos='http://cran.us.r-project.org')
#install.packages('gridExtra', repos='http://cran.us.r-project.org')
#quit(save = "default")

# Make enrichment plots
sbatch $CODE/enrichment/launch_enrichmentplot.sh $CODE $FIRST_SAMPLE $LAST_SAMPLE 5 $DATA/counts/CX${CTRL_SAMPLE}_n5_mapq3
0_sorted_qual30.count $CHROM_NAME $DATA/enrichment $DATA/counts
