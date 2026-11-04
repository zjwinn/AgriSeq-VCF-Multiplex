#!/bin/bash
#SBATCH --job-name=AgriSeq_Subset_Sort
#SBATCH --output=agriseq_subset_sort_%j.out
#SBATCH --error=agriseq_subset_sort_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=medium
#SBATCH --account=your_project_account
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.name@usda.gov

# Print start time and job details
echo "Job started on $(hostname) at $(date)"

# Determine repository root directory
# (When run via sbatch, SLURM_SUBMIT_DIR preserves the directory where sbatch was called)
if [ -n "$SLURM_SUBMIT_DIR" ]; then
    REPO_DIR="$SLURM_SUBMIT_DIR"
else
    REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

# Environment Configuration (Choose Method 1 or Method 2):
# Method 1 (Recommended on Ceres): Use pre-installed cluster modules
module purge
module load bcftools
module load htslib

# Method 2 (Alternative): Activate custom Conda environment on Ceres
# module load miniconda3
# conda activate agriseq_vcf_env

# Define input parameters (should be realpath of each input)
INPUT_VCF="Merged_Dataset.vcf.gz"
OUTPUT_VCF="Subset_Filtered_Sorted.vcf.gz"
TARGETS_FILE="$REPO_DIR/target_manifests/Wheat_AgSeqv2_Targets_Manifest.txt" # Or example: $REPO_DIR/example_inputs/example_target_manifest.txt
LINES_FILE="$REPO_DIR/example_inputs/example_lines_list.txt"
SORT_MEMORY="28g"
BIALLELIC_FLAG="-b" # Set to "-b" for strict biallelic REF/ALT compliance (T3/Breedbase, GWAS, GS), or "" for non-destructive mode

# Run the subsetting and sorting pipeline
"$REPO_DIR/AgriSeq_Subset_and_Sort.sh" \
  --input "$INPUT_VCF" \
  --output "$OUTPUT_VCF" \
  --targets "$TARGETS_FILE" \
  --lines "$LINES_FILE" \
  --memory "$SORT_MEMORY" \
  $BIALLELIC_FLAG

# Print finish time
echo "Job completed at $(date)"
