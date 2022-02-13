# SNPsplitPipeline
To search for parent-of-origin effects on gene expression, our lab has studied hybrid mice. To study parent-of-origin effects on chromotin confirmation and gene-regulationn, we have studied these same hybrid mice with PLAC-seq--a HiC derivitive. The variants that distinguish the maternally and paternally inheritied DNA alow us to identify contacts that differ with mouse strain and parent-of-origin effects. This pipeline is built to align hybrid mouse PLAC-seq data and sort the reads that contain varaints specific to each parent genome. These reads are then aligned seperatly by the MAPS pipeline. The MAPS pipeline calculates the significance of contacts for both alleles.

I run this pipeline using the SLURM scheduler and make use of 'modules' to load programs including trimmomatic, bowtie2, samtools and bedtools. The script should be modified for different systems.

# Prerequisites
samtools
bedtools
bowtie2
trimmomatic
MAPS

Install the MAPS pipeline locally. See the MAPS page for the required dependancies (https://github.com/ijuric/MAPS).
Also install SNPsplit locally (https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/). SNPsplit requires a locally samtools installation.

These data come from PLAC-seq experiments run with tissue from C57/Cast hybrid mice bred to allow maternal and paternal RNA-seq reads to be distinguished after alignment. Trim and align PLAC-seq FASTQ files. I've set the parameters for 150 bp paired-end reads.
Use known C57/Cast SNPS to sort reads into separate Cast and C57 FASTQ files. Use a bowtie alignment index wiht N masked Variants. Pass the separate FASTQ files to the MAPS pipeline to call significant cis contacts for each allele seperatly. MAPS parameters can be changed in the shell script 'run_pipeline_cl.sh'. Move run_pipeline_cl.sh to ~/MAPS/bin/. This will ultimately give us the data we need to visualize loops for the two alleles.
Requres a variant N-masked bowtie2 index file for alignment. Use bowtie2-build
The MAPS pipeline requires of ChIP-seq peaks BED file from analgous samples.
