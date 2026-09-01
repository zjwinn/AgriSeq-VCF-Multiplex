#!/bin/bash

# ==============================================================================
# Example Wrapper Script for Subsetting and Sorting
# Demonstrates calling AgriSeq_Subset_and_Sort.sh with command line flags.
# ==============================================================================

# Determine repository root directory
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# 1. Activate the required Conda environment
echo "Activating conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || eval "$(conda shell.bash hook 2>/dev/null)"
conda activate agriseq_vcf_env

# 2. Define input parameters using public example templates
INPUT_VCF="Merged_Dataset.vcf.gz"
OUTPUT_VCF="Subset_Filtered_Sorted.vcf.gz"
TARGETS_FILE="$REPO_DIR/example_target_manifest.txt"
LINES_FILE="$REPO_DIR/example_lines_list.txt"
SORT_MEMORY="16g"

# 3. Execute the subsetting and sorting pipeline (again check input names to make sure they are not just generic inputs)
echo "Executing AgriSeq_Subset_and_Sort.sh pipeline... (Output is being saved to subset_and_sort.log)"
"$REPO_DIR/AgriSeq_Subset_and_Sort.sh" \
  --input "$INPUT_VCF" \
  --output "$OUTPUT_VCF" \
  --targets "$TARGETS_FILE" \
  --lines "$LINES_FILE" \
  --memory "$SORT_MEMORY" &> subset_and_sort.log

echo "Job completed. Log written to subset_and_sort.log"
