#!/bin/bash
"""

Scripts and source data to reproduce figures of the manuscript:
Hess, Anja et al. 2025: Non-disruptive monitoring of cellular dynamics with cell-free DNA methylation

Max Planck Institute for Molecular Genetics, Berlin, Germany

Author: Anja Hess
Date: 2025 APRIL 07

"""

####################################################################################
# BEFORE YOU START
####################################################################################
# Step 1
# For bam density plots please download the fastq.gz files from GEO (GSE293866),
# align to mm10 (details under Data processing on GSE293866).
# If you have slightly different alignment #
# pipelines this should not drastically affect read density. E
# ach bam should be around about ~50GB.

# Step 2
# Put them into this folder. You may adjust names for BAM1 and BAM2 below.

####################################################################################
# FIG 1E - Karyplot for read density
####################################################################################
# Bash script to generate coverage bigwigs for all bam files in directory
# USGAE: bash Fig_1H_karyo.sh
# AUTHOR: Anja Hess
# DATE: 20241212

R_script_loc="./utils/karyo.R"

# Attention! Names may differ depd. on your pipeline
BAM1="mpimg_L31945-1_AM-WGBS-622_S8_mm10.bsmap.srt.rd.bam"
BAM2="mpimg_L31953-1_AM-WGBS-630_S16_mm10.bsmap.srt.rd.bam"

# 19 autosomes.. call the R script for each
for i in {1..19}; do
  (echo $FILE
  Rscript $R_script_loc $BAM1 $BAM2 "chr${i}" mm10 "chr${i}_mm10.pdf")&
done
wait


# END OF SCRIPT