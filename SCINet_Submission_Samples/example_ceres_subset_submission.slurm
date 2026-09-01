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
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Source Conda from the default user Miniconda installation
source ~/miniconda3/etc/profile.d/conda.sh
conda activate agriseq_vcf_env

# Define input parameters (should be realpath of each input)
INPUT_VCF="Merged_Dataset.vcf.gz"
OUTPUT_VCF="Subset_Filtered_Sorted.vcf.gz"
TARGETS_FILE="$REPO_DIR/example_target_manifest.txt"
LINES_FILE="$REPO_DIR/example_lines_list.txt"
SORT_MEMORY="28g"

# Run the subsetting and sorting pipeline
"$REPO_DIR/AgriSeq_Subset_and_Sort.sh" \
  --input "$INPUT_VCF" \
  --output "$OUTPUT_VCF" \
  --targets "$TARGETS_FILE" \
  --lines "$LINES_FILE" \
  --memory "$SORT_MEMORY"

# Print finish time
echo "Job completed at $(date)"
