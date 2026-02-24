#!/bin/bash
#===============================================================================
# run_pipeline_cl.sh
#
# Wrapper script for the MAPS pipeline (MAPS step only, feather disabled).
#
# Description:
#   Runs the MAPS (Model-based Analysis of PLAC-seq) pipeline on allele-split
#   FASTQ files produced by SNPsplitPipeline.sh. This script handles binning,
#   regression, and peak calling to identify significant cis chromatin contacts.
#
# Prerequisites:
#   - Python (with pysam, pybedtools installed)
#   - R / Rscript
#   - samtools and bedtools in $PATH
#   - MAPS pipeline (https://github.com/HuMingLab/MAPS)
#   - MAPS_data_files for the target organism (genomic features, BWA index)
#
# Usage:
#   run_pipeline_cl.sh <dataset_name> <fastq_dir> <outdir> <macs2_filepath>
#
# Arguments:
#   dataset_name     Sample/dataset identifier
#   fastq_dir        Directory containing input FASTQ files
#                    (expects <dataset_name>_R1.fastq and _R2.fastq)
#   outdir           Output directory for MAPS results
#   macs2_filepath   Path to ChIP-seq peaks BED file (from MACS2)
#
# Notes:
#   - Copy this script to ~/MAPS/bin/ before running SNPsplitPipeline.sh
#   - Adjust MAPS parameters in the configuration section below
#
# Author:  Elliott Ferris
# Date:    2020-08-10
#===============================================================================

set -euo pipefail

#---------------------------------------
# Parse command-line arguments
#---------------------------------------
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <dataset_name> <fastq_dir> <outdir> <macs2_filepath>"
    exit 1
fi

dataset_name="$1"
fastq_dir="$2"
outdir="$3"
macs2_filepath="$4"

#===============================================================================
# MAPS Configuration
# ** Edit these parameters as needed **
#===============================================================================

# Tool paths
python_path=~/anaconda3/bin/python
Rscript_path="/usr/bin/Rscript"       # <-- Set your Rscript path

# Pipeline control
feather=1                              # 1 = run feather step; 0 = skip
maps=1                                 # 1 = run MAPS step;   0 = skip
number_of_datasets=1

# Organism and reference genome
organism="mm10"                        # Options: mm9, mm10, hg19, hg38
bwa_index=""                           # <-- Set your BWA index path (leave empty for default)

# MAPS parameters
fastq_format=".fastq"
bin_size=10000
binning_range=1000000
fdr=2                                  # Used for labeling only; do not change
filter_file="None"
generate_hic=1
mapq=30
length_cutoff=1000
threads=54
model="pospoisson"                     # Options: "pospoisson", "negbinom"
sex_chroms_to_process="X"              # Options: "X", "Y", "XY", or anything else for autosomes only
optical_duplicate_distance=0

#===============================================================================
# Multi-dataset merging (only used when number_of_datasets > 1)
#===============================================================================
dataset1=""
dataset2=""
dataset3=""
dataset4=""

# Override feather output symlink (only if feather=0 and using specific output)
feather_output_symlink=""

#===============================================================================
# Derived variables — do not modify
#===============================================================================

DATE=$(date '+%Y%m%d_%H%M%S')
cwd="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null && pwd)"

fastq1="${fastq_dir}/${dataset_name}_R1${fastq_format}"
fastq2="${fastq_dir}/${dataset_name}_R2${fastq_format}"

feather_output="${outdir}/feather_output/${dataset_name}_${DATE}"
if [ -z "${feather_output_symlink}" ]; then
    feather_output_symlink="${outdir}/feather_output/${dataset_name}_current"
fi

resolution=$((bin_size / 1000))
per_chr="True"  # Set to "False" to skip per-chromosome BED/BEDPE output
feather_logfile="${feather_output}/${dataset_name}.feather.log"
hic_dir="tempfiles/hic_tempfiles"

#---------------------------------------
# Set organism-specific paths
#---------------------------------------
chr_count=0
case "${organism}" in
    mm10)
        [ -z "${bwa_index}" ] && bwa_index="${cwd}/../MAPS_data_files/${organism}/BWA_index/mm10_chrAll.fa"
        genomic_feat_filepath="${cwd}/../MAPS_data_files/${organism}/genomic_features/F_GC_M_MboI_${resolution}Kb_el.mm10.txt"
        chr_count=19
        ;;
    mm9)
        [ -z "${bwa_index}" ] && bwa_index="${cwd}/../MAPS_data_files/${organism}/BWA_index/mm9.fa"
        genomic_feat_filepath="${cwd}/../MAPS_data_files/${organism}/genomic_features/F_GC_M_MboI_${resolution}Kb_el.mm9.txt"
        chr_count=19
        ;;
    hg19)
        [ -z "${bwa_index}" ] && bwa_index="${cwd}/../MAPS_data_files/${organism}/BWA_index/hg19.fa"
        genomic_feat_filepath="${cwd}/../MAPS_data_files/${organism}/genomic_features/F_GC_M_MboI_${resolution}Kb_el.hg19.txt"
        chr_count=22
        ;;
    hg38)
        [ -z "${bwa_index}" ] && bwa_index="${cwd}/../MAPS_data_files/${organism}/BWA_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
        genomic_feat_filepath="${cwd}/../MAPS_data_files/${organism}/genomic_features/F_GC_M_MboI_${resolution}Kb_el.GRCh38.txt"
        chr_count=22
        ;;
    *)
        echo "ERROR: Unsupported organism '${organism}'. Use mm9, mm10, hg19, or hg38."
        exit 1
        ;;
esac

#---------------------------------------
# Process sex chromosome settings
#---------------------------------------
if [[ "${sex_chroms_to_process}" != "X" && \
      "${sex_chroms_to_process}" != "Y" && \
      "${sex_chroms_to_process}" != "XY" ]]; then
    sex_chroms_to_process="NA"
    sex_chroms=""
else
    sex_chroms="${sex_chroms_to_process}"
fi

long_bedpe_dir="${feather_output_symlink}/"
short_bed_dir="${feather_output_symlink}/"
maps_output="${outdir}/MAPS_output/${dataset_name}_${DATE}/"
maps_output_symlink="${outdir}/MAPS_output/${dataset_name}_current"

#===============================================================================
# Run MAPS
#===============================================================================

if [ "${maps}" -eq 1 ]; then
    echo "============================================"
    echo " MAPS Pipeline"
    echo " Dataset: ${dataset_name}"
    echo " Output:  ${maps_output}"
    echo " Started: $(date)"
    echo "============================================"

    mkdir -p "${maps_output}"

    echo "[$(date)] Generating MAPS runfile..."
    ${python_path} "${cwd}/MAPS/make_maps_runfile.py" \
        "${dataset_name}" \
        "${maps_output}" \
        "${macs2_filepath}" \
        "${genomic_feat_filepath}" \
        "${long_bedpe_dir}" \
        "${short_bed_dir}" \
        "${bin_size}" \
        "${chr_count}" \
        "${maps_output}" \
        "${sex_chroms_to_process}" \
        --BINNING_RANGE "${binning_range}"

    echo "[$(date)] Running MAPS..."
    ${python_path} "${cwd}/MAPS/MAPS.py" \
        "${maps_output}maps_${dataset_name}.maps"

    echo "[$(date)] Running regression and peak calling..."
    ${Rscript_path} "${cwd}/MAPS/MAPS_regression_and_peak_caller.r" \
        "${maps_output}" \
        "${dataset_name}.${resolution}k" \
        "${bin_size}" \
        "${chr_count}${sex_chroms}" \
        "${filter_file}" \
        "${model}"

    echo "[$(date)] Formatting peaks..."
    ${Rscript_path} "${cwd}/MAPS/MAPS_peak_formatting.r" \
        "${maps_output}" \
        "${dataset_name}.${resolution}k" \
        "${fdr}" \
        "${bin_size}"

    # Archive a copy of this script with the output
    cp "$(readlink -f "$0")" "${maps_output}/execution_script_copy"
    chmod 755 "${maps_output}"
    ln -sfn "${maps_output}" "${maps_output_symlink}"

    echo "============================================"
    echo " MAPS complete: $(date)"
    echo "============================================"
fi
