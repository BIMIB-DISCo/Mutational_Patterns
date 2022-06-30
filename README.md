# Mutational_Patterns

## Installation
Install the following R packages to run the R script:
* MutationalPatterns
* BSgenome
* NMF
* Biostrings
* stringr
* stringdist
* Matching
* gplots

## Data
From the folder where the R script file is located, create the "vcf_organoids" folder.
Download all vcf_organoids.7z files and decompress them in a "vcf_organoids" folder.

## Script
Execute the R script

## Data preparation
If a vcf file contains variant calling from multiple samples and 
the header does not contains the contigs for the chromosomes 
(as some of those in the accelerator repository),
then split_vcf.sh is the right script for you.
It splits vcf per sample and it makes sure chr contigs are in the header.
Contigs size were taken from parental organoid vcf file (2019).
This is another place where to find contigs: https://www.ncbi.nlm.nih.gov/grc/human/data
