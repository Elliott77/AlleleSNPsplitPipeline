# SNPsplitPipeline
I run this pipeline on the University of Utah's CHPC system using the SLURM scheduler and make use of 'module' to load programs including trimmomatic, bowtie2, samtools and bedtools.

# Prerequisites
samtools

bedtools

bowtie2

trimmomatic.

Install the software MAPS locally. See the MAPS page for the required dependancies (https://github.com/ijuric/MAPS).
I also have the program SNPsplit installed locally (https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/). SNPsplit requires a locally samtools installation.

These data come from PLAC-seq experiments run with tissue from C57/Cast hybrid mice bred to allow maternally and paternally inherited reads to be distinguished. Trim and align PLAC-seq FASTQ files. I've set the parameters for 150 bp paired-end reads.
Use known C57/Cast SNPS to sort reads into separate Cast and C57 FASTQ files. Use a bowtie alignment index wiht N masked CAST Variants. Pass the separate FASTQ to the MAPS pipeline to call significant cis contacts for each allele seperatly. MAPS parameters can be changed in the shell script 'run_pipeline_cl.sh'. This will ultimately give us the data we need to visualize loops that differ between the two alleles.
Requres a Cast variant N-masked bowtie2 index file for alignment
MAPS requires of ChIP-seq peaks BED file from analgous samples.
