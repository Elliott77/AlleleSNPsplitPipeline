#!/bin/sh
## SNPsplitPipeline.sh
## Elliott Ferris
## 8/10/2020

## Useage: SNPsplitPipeline.sh <fastq_prefix> <fastq_directory> <Adaptor sequence>  <threads>

#SBATCH --time=18:50:00 # Walltime
#SBATCH --nodes=1          
#SBATCH --ntasks=1       
#SBATCH --account=gregg
#SBATCH --partition=<partion-name>
#SBATCH --account=<account>
#SBATCH --job-name=SNP_split_$1
#SBATCH -o <log files path>/slurm_snpsplit-%j.out-%N
#SBATCH -e <log files path>/stderr-%j.txt
## useage
## ./SNPsplitPiplineTrim2 FastQprefix /FastQDir ACGTCGTC

FASTQ_PREFIX=$1
threads=$4
fastq_dir=$2
scratch_sam=<scratch directory for SAM files>
echo $threads
ADAPTOR=$3
AlleleFASQ_Directory=<<allele fastq directory>> 
SNPFileDir=<SNP FIle Directory>
GRCm38_68_cast_N_masked=<GRCm38 cast N masked alignment index>

cd ${fastq_dir}
echo ">PrefixPE/1" > ${FASTQ_PREFIX}_PE.fa
## Illumina adaptor
echo AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT >> ${FASTQ_PREFIX}_PE.fa
echo ">PrefixPE/2" >> ${FASTQ_PREFIX}_PE.fa
echo ${ADAPTOR} >> ${FASTQ_PREFIX}_PE.fa

## trim reads
module load trimmomatic
cd ${fastq_dir}
trimmomatic PE -threads $threads ${FASTQ_PREFIX}_R1.fastq.gz ${FASTQ_PREFIX}_R2.fastq.gz ${FASTQ_PREFIX}_paired_R1.fastq ${FASTQ_PREFIX}_unpaired_R1.fastq ${FASTQ_PREFIX}_paired_R2.fastq ${FASTQ_PREFIX}_unpaired_R2.fastq ILLUMI
NACLIP:{FASTQ_PREFIX}_PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
cd ${fastq_dir}

module unload trimmomatic
echo "Trimmomatic done at `date`"##
rm  ${FASTQ_PREFIX}_unpaired_R1.fastq  ${FASTQ_PREFIX}_unpaired_R2.fastq.gz

macs2=<ChIP-seq peaks BED file>

## align reads with bowtie
module load bowtie2
echo "align reads at `date`"
bowtie2 -q --sensitive --threads $threads  -x $GRCm38_68_cast_N_masked --end-to-end -1 ${fastq_dir}/${FASTQ_PREFIX}_paired_R1.
fastq -2 ${fastq_dir}/${FASTQ_PREFIX}_paired_R2.fastq -S ${scratch_sam}${FASTQ_PREFIX}.sam

module unload bowtie2
rm ${FASTQ_PREFIX}_paired_R1.fastq.gz ${FASTQ_PREFIX}_paired_R2.fastq.gz
## Sort SAM to BAM
echo "Sort SAM to BAM at `date`"

module load samtools
STP=$(which samtools)
cd $scratch_sam
samtools sort -@ $threads -O BAM -o ${FASTQ_PREFIX}.bam ${FASTQ_PREFIX}.sam
rm {FASTQ_PREFIX}.sam
## Sort Reads
echo "Starting SNPsplit at `date`"
cd $scratch_sam
## use the program 'SNPsplit' to sort reads as containing C57 or Cast specific variants.
SNPsplit --samtools_path $STP --conflicting --paired --snp_file ${SNPFileDir}/all_SNPs_CAST_EiJ_GRCm38.txt.gz ${FASTQ_PREFIX}.bam
echo "SNPsplit done at `date`"

## Convert BAM to FASTQ
echo "BAM to FASTQ"

module load bedtools

cd $scratch_sam
## paired end FASTQ reads
samtools view -hbf 64 ${FASTQ_PREFIX}.genome1.bam | bedtools bamtofastq -i stdin -fq ${AlleleFASQ_Directory}/c57/${FASTQ_PREFIX}_c57_R1.fastq
samtools view -hbf 128 ${FASTQ_PREFIX}.genome1.bam | bedtools bamtofastq -i stdin -fq ${AlleleFASQ_Directory}/c57/${FASTQ_PREFIX}_c57_R2.fastq
samtools view -hbf 64 ${FASTQ_PREFIX}.genome2.bam | bedtools bamtofastq -i stdin -fq ${AlleleFASQ_Directory}/cast/${FASTQ_PREFIX}_cast_R1.fastq
samtools view -hbf 128 ${FASTQ_PREFIX}.genome2.bam | bedtools bamtofastq -i stdin -fq ${AlleleFASQ_Directory}/cast/${FASTQ_PREFIX}_cast_R2.fastq

module unload bedtools samtools
echo "BAM to FASTQ Finished, Time for MAPS"
echo "Starting MAPS at `date`"
module load R/3.6.1
module load bedtools samtools/1.10 bwa
cd ${MAPS_Output}
mkdir ${FASTQ_PREFIX}_c57
mkdir ${FASTQ_PREFIX}_cast
cd ~/MAPS/bin
./run_pipeline_cl.sh ${FASTQ_PREFIX}_c57 ${AlleleFASQ_Directory}/c57 ${MAPS_Output}/${FASTQ_PREFIX}_c57 ${macs2}
./run_pipeline_cl.sh ${FASTQ_PREFIX}_cast ${AlleleFASQ_Directory}/cast ${MAPS_Output}/${FASTQ_PREFIX}_cast  ${macs2}
echo "MAPS Done!"
cd $scratch_sam
rm ${FASTQ_PREFIX}.sam
rm ${FASTQ_PREFIX}.genome1.bam 
rm ${FASTQ_PREFIX}.genome2.bam 
cd ${AlleleFASQ_Directory}
rm ${FASTQ_PREFIX}_c57_R1.fastq ${FASTQ_PREFIX}_c57_R2.fastq ${FASTQ_PREFIX}_cast_R1.fastq ${FASTQ_PREFIX}_cast_R2.fastq

module purge
echo "End of program at `date`")

