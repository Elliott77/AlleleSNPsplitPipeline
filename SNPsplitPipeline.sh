#!/bin/bash
#===============================================================================
# SNPsplitPipeline.sh
#
# Allele-specific Hi-C analysis pipeline using SNPsplit and MAPS.
#
# Description:
#   Uses known SNPs to sort paired-end reads into separate Strain A (C57)
#   and Strain B (CAST) FASTQ files via SNPsplit. Reads are aligned to a
#   variant N-masked bowtie2 index. The allele-split FASTQs are then passed
#   to the MAPS pipeline to call significant cis contacts for each allele
#   separately, yielding data for allele-specific loop visualization.
#
# Prerequisites:
#   - trimmomatic
#   - bowtie2  (+ a variant N-masked bowtie2 index; build with bowtie2-build)
#   - samtools
#   - bedtools
#   - SNPsplit  (Babraham Bioinformatics)
#   - MAPS pipeline (https://github.com/HuMingLab/MAPS)
#     * Copy run_pipeline_cl.sh to ~/MAPS/bin/ and configure MAPS parameters
#       within that script.
#   - A ChIP-seq peaks BED file from analogous samples (for MAPS)
#   - SNP file for SNPsplit (e.g., all_SNPs_CAST_EiJ_GRCm38.txt.gz)
#
# Usage:
#   SNPsplitPipeline.sh <fastq_prefix> <fastq_directory> <adaptor_sequence> <threads>
#
# Arguments:
#   fastq_prefix       Prefix for input FASTQ files (expects <prefix>_R1.fastq.gz
#                      and <prefix>_R2.fastq.gz in the FASTQ directory)
#   fastq_directory    Path to directory containing input FASTQ files
#   adaptor_sequence   3' adaptor sequence for trimming
#   threads            Number of threads for parallel processing
#
# Author:  Elliott Ferris
# Date:    2020-08-10
#===============================================================================

set -euo pipefail

#---------------------------------------
# Parse command-line arguments
#---------------------------------------
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <fastq_prefix> <fastq_directory> <adaptor_sequence> <threads>"
    exit 1
fi

FASTQ_PREFIX="$1"
FASTQ_DIR="$2"
ADAPTOR="$3"
THREADS="$4"

#---------------------------------------
# User-configurable paths
# ** Edit these before running **
#---------------------------------------

# Scratch directory for intermediate SAM/BAM files
SCRATCH_SAM="/path/to/scratch/sam/"

# Output directory for allele-split FASTQ files (must contain c57/ and cast/ subdirs)
ALLELE_FASTQ_DIR="/path/to/allele_fastq/"

# SNPsplit SNP file (e.g., all_SNPs_CAST_EiJ_GRCm38.txt.gz)
SNP_FILE="/path/to/snp_files/all_SNPs_CAST_EiJ_GRCm38.txt.gz"

# Bowtie2 index prefix for the variant N-masked reference genome
BOWTIE2_INDEX="/path/to/GRCm38_68_cast_N_masked"

# ChIP-seq peaks BED file (required by MAPS)
CHIPSEQ_PEAKS="/path/to/chipseq_peaks.bed"

# MAPS output directory
MAPS_OUTPUT="/path/to/maps_output/"

#===============================================================================
# Pipeline
#===============================================================================

echo "============================================"
echo " SNPsplit Pipeline"
echo " Sample:  ${FASTQ_PREFIX}"
echo " Threads: ${THREADS}"
echo " Started: $(date)"
echo "============================================"

#---------------------------------------
# 1. Build adaptor FASTA for Trimmomatic
#---------------------------------------
echo "[$(date)] Building adaptor FASTA..."
cd "${FASTQ_DIR}"

cat > "${FASTQ_PREFIX}_PE.fa" <<EOF
>PrefixPE/1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
${ADAPTOR}
EOF

#---------------------------------------
# 2. Trim reads with Trimmomatic
#---------------------------------------
echo "[$(date)] Running Trimmomatic..."

trimmomatic PE \
    -threads "${THREADS}" \
    "${FASTQ_PREFIX}_R1.fastq.gz" \
    "${FASTQ_PREFIX}_R2.fastq.gz" \
    "${FASTQ_PREFIX}_paired_R1.fastq" \
    "${FASTQ_PREFIX}_unpaired_R1.fastq" \
    "${FASTQ_PREFIX}_paired_R2.fastq" \
    "${FASTQ_PREFIX}_unpaired_R2.fastq" \
    "ILLUMINACLIP:${FASTQ_PREFIX}_PE.fa:2:30:10:2:keepBothReads" \
    LEADING:3 \
    TRAILING:3 \
    MINLEN:36

echo "[$(date)] Trimmomatic complete."

# Clean up unpaired reads and adaptor file
rm -f "${FASTQ_PREFIX}_unpaired_R1.fastq" \
      "${FASTQ_PREFIX}_unpaired_R2.fastq" \
      "${FASTQ_PREFIX}_PE.fa"

#---------------------------------------
# 3. Align reads with Bowtie2
#---------------------------------------
echo "[$(date)] Aligning reads with Bowtie2..."

bowtie2 -q \
    --sensitive \
    --end-to-end \
    --threads "${THREADS}" \
    -x "${BOWTIE2_INDEX}" \
    -1 "${FASTQ_DIR}/${FASTQ_PREFIX}_paired_R1.fastq" \
    -2 "${FASTQ_DIR}/${FASTQ_PREFIX}_paired_R2.fastq" \
    -S "${SCRATCH_SAM}/${FASTQ_PREFIX}.sam"

echo "[$(date)] Alignment complete."

# Clean up trimmed paired reads
rm -f "${FASTQ_PREFIX}_paired_R1.fastq" \
      "${FASTQ_PREFIX}_paired_R2.fastq"

#---------------------------------------
# 4. Sort SAM to BAM
#---------------------------------------
echo "[$(date)] Sorting SAM to BAM..."

cd "${SCRATCH_SAM}"
samtools sort \
    -@ "${THREADS}" \
    -O BAM \
    -o "${FASTQ_PREFIX}.bam" \
    "${FASTQ_PREFIX}.sam"

rm -f "${FASTQ_PREFIX}.sam"

#---------------------------------------
# 5. Split reads by allele with SNPsplit
#---------------------------------------
echo "[$(date)] Running SNPsplit..."

STP=$(which samtools)

SNPsplit \
    --samtools_path "${STP}" \
    --conflicting \
    --paired \
    --snp_file "${SNP_FILE}" \
    "${FASTQ_PREFIX}.bam"

echo "[$(date)] SNPsplit complete."

#---------------------------------------
# 6. Convert allele-split BAMs to FASTQ
#---------------------------------------
echo "[$(date)] Converting BAM to FASTQ..."

mkdir -p "${ALLELE_FASTQ_DIR}/c57" "${ALLELE_FASTQ_DIR}/cast"

# C57 (genome1) — Read 1 and Read 2
samtools view -hbf 64  "${FASTQ_PREFIX}.genome1.bam" | \
    bedtools bamtofastq -i stdin -fq "${ALLELE_FASTQ_DIR}/c57/${FASTQ_PREFIX}_c57_R1.fastq"

samtools view -hbf 128 "${FASTQ_PREFIX}.genome1.bam" | \
    bedtools bamtofastq -i stdin -fq "${ALLELE_FASTQ_DIR}/c57/${FASTQ_PREFIX}_c57_R2.fastq"

# CAST (genome2) — Read 1 and Read 2
samtools view -hbf 64  "${FASTQ_PREFIX}.genome2.bam" | \
    bedtools bamtofastq -i stdin -fq "${ALLELE_FASTQ_DIR}/cast/${FASTQ_PREFIX}_cast_R1.fastq"

samtools view -hbf 128 "${FASTQ_PREFIX}.genome2.bam" | \
    bedtools bamtofastq -i stdin -fq "${ALLELE_FASTQ_DIR}/cast/${FASTQ_PREFIX}_cast_R2.fastq"

echo "[$(date)] BAM to FASTQ conversion complete."

#---------------------------------------
# 7. Run MAPS for each allele
#---------------------------------------
echo "[$(date)] Starting MAPS pipeline..."

mkdir -p "${MAPS_OUTPUT}/${FASTQ_PREFIX}_c57"
mkdir -p "${MAPS_OUTPUT}/${FASTQ_PREFIX}_cast"

cd ~/MAPS/bin

./run_pipeline_cl.sh \
    "${FASTQ_PREFIX}_c57" \
    "${ALLELE_FASTQ_DIR}/c57" \
    "${MAPS_OUTPUT}/${FASTQ_PREFIX}_c57" \
    "${CHIPSEQ_PEAKS}"

./run_pipeline_cl.sh \
    "${FASTQ_PREFIX}_cast" \
    "${ALLELE_FASTQ_DIR}/cast" \
    "${MAPS_OUTPUT}/${FASTQ_PREFIX}_cast" \
    "${CHIPSEQ_PEAKS}"

echo "[$(date)] MAPS complete."

#---------------------------------------
# 8. Clean up intermediate files
#---------------------------------------
echo "[$(date)] Cleaning up..."

rm -f "${SCRATCH_SAM}/${FASTQ_PREFIX}.bam"
rm -f "${SCRATCH_SAM}/${FASTQ_PREFIX}.genome1.bam"
rm -f "${SCRATCH_SAM}/${FASTQ_PREFIX}.genome2.bam"

rm -f "${ALLELE_FASTQ_DIR}/c57/${FASTQ_PREFIX}_c57_R1.fastq" \
      "${ALLELE_FASTQ_DIR}/c57/${FASTQ_PREFIX}_c57_R2.fastq" \
      "${ALLELE_FASTQ_DIR}/cast/${FASTQ_PREFIX}_cast_R1.fastq" \
      "${ALLELE_FASTQ_DIR}/cast/${FASTQ_PREFIX}_cast_R2.fastq"

echo "============================================"
echo " Pipeline complete: $(date)"
echo "============================================"

