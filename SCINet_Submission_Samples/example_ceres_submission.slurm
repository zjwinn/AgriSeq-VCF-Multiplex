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
#SBATCH --account=your_project_account
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.name@usda.gov

# Print start time and job details
echo "Job started on $(hostname) at $(date)"

# Determine repository root directory
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Source Conda from the default user Miniconda installation
source ~/miniconda3/etc/profile.d/conda.sh
conda activate agriseq_vcf_env

# Define input parameters (should be realpath of any input files)
ZIP_DIR="raw_zip_files"
REF_INDEX="CHS21_OT"
CONTIGS="Chr1A,Chr1B,Chr1D,Chr2A,Chr2B,Chr2D,Chr3A,Chr3B,Chr3D,Chr4A,Chr4B,Chr4D,Chr5A,Chr5B,Chr5D,Chr6A,Chr6B,Chr6D,Chr7A,Chr7B,Chr7D"
OUTPUT_VCF="Merged_Dataset.vcf.gz" # this is a generic location and name 
MAP_FILE="$REPO_DIR/example_sample_map.txt" # This is an example file replace with a real one
BATCH_SIZE=100
NON_NATIVE_PREFIX="OT"
DEDUP_FLAG="-d" # This is an option to resolve marker names, usually leave this alone

# Run the AgriSeq merger pipeline
"$REPO_DIR/AgriSeq_Merger.sh" \
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
