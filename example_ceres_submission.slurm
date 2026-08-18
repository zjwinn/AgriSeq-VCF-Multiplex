#!/bin/bash
#SBATCH --job-name=AgriSeq_Merger
#SBATCH --output=agriseq_merge_%j.out
#SBATCH --error=agriseq_merge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=medium
#SBATCH --account=breeding_hwwgru
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.name@usda.gov

# Print start time and job details
echo "Job started on $(hostname) at $(date)"

# Source Conda from the default user Miniconda installation
source ~/miniconda3/etc/profile.d/conda.sh
conda activate agriseq_vcf_env

# Define input parameters
ZIP_DIR="NAG_Zip_Files"
REF_INDEX="CHS21_OT"
CONTIGS="Chr1A,Chr1B,Chr1D,Chr2A,Chr2B,Chr2D,Chr3A,Chr3B,Chr3D,Chr4A,Chr4B,Chr4D,Chr5A,Chr5B,Chr5D,Chr6A,Chr6B,Chr6D,Chr7A,Chr7B,Chr7D"
OUTPUT_VCF="USDA_AgriSeq_2025_2026.vcf.gz"
MAP_FILE="KSM2026-GYT_Sample_Names.txt"
BATCH_SIZE=100
NON_NATIVE_PREFIX="OT"
DEDUP_FLAG="-d"

# Run the AgriSeq merger pipeline
./AgriSeq_Merger.sh \
  --zip-dir "$ZIP_DIR" \
  --ref-index "$REF_INDEX" \
  --contigs "$CONTIGS" \
  --output "$OUTPUT_VCF" \
  --batch-size "$BATCH_SIZE" \
  --map-file "$MAP_FILE" \
  --non-native "$NON_NATIVE_PREFIX" \
  $DEDUP_FLAG

# Print finish time
echo "Job completed at $(date)"
