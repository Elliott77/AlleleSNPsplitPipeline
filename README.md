# SNPsplitPipeline
To search for parent-of-origin and strain effects on gene expression, our lab has studied hybrid mice. To study parent-of-origin and strain effects on chromotin confirmation and gene-regulationn, we have studied these same hybrid mice with PLAC-seq--a HiC derivitive. The variants that distinguish the maternally and paternally inheritied DNA alow us to identify contacts that differ with mouse strain and parent-of-origin effects. This pipeline is built to align hybrid mouse PLAC-seq data and sort the reads that contain varaints specific to each parent genome. These reads are then aligned seperatly by the MAPS pipeline. The MAPS pipeline calculates the significance of contacts for both alleles.

# Prerequisites
samtools
bedtools
bowtie2
trimmomatic
MAPS

Install the MAPS pipeline locally. See the MAPS page for the required dependancies (https://github.com/ijuric/MAPS).
Also install SNPsplit locally (https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/). SNPsplit requires a locally samtools installation.

These data come from PLAC-seq experiments run with tissue from C57/Cast hybrid mice bred to allow maternal and paternal RNA-seq reads to be distinguished after alignment. Trim and align PLAC-seq FASTQ files. I've set the parameters for 150 bp paired-end reads. Modify the script as needed.
Use known SNPS to sort reads into separate Strain A and Strain B FASTQ files. Use a bowtie alignment index wiht N masked Variants. Pass the separate FASTQ files to the MAPS pipeline to call significant cis contacts for each allele seperatly. MAPS parameters can be changed in the script 'run_pipeline_cl.sh'. Move run_pipeline_cl.sh to ~/MAPS/bin/. This will ultimately yield the data needed to visualize loops for the two alleles.
Requres a variant N-masked bowtie2 index file for alignment. Use bowtie2-build to create
The MAPS pipeline requires of ChIP-seq peaks BED file from analgous samples.
