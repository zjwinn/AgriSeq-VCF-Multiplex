#!/bin/bash

# ==============================================================================
# Example Wrapper Script for AgriSeq Variant Merging
# Demonstrates calling AgriSeq_Merger.sh with command line flags.
# ==============================================================================

# Determine repository root directory
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# 1. Activate the required Conda environment
echo "Activating conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || eval "$(conda shell.bash hook 2>/dev/null)"
conda activate agriseq_vcf_env

# 2. Execute the pipeline using generic example parameters (again check you inputs to make sure they are correct and not examples)
echo "Executing AgriSeq_Merger.sh pipeline... (Output is being saved to pipeline.log)"
"$REPO_DIR/AgriSeq_Merger.sh" \
  --zip-dir "raw_zip_files" \
  --ref-index "CHS21_OT" \
  --contigs "Chr1A,Chr1B,Chr1D,Chr2A,Chr2B,Chr2D,Chr3A,Chr3B,Chr3D,Chr4A,Chr4B,Chr4D,Chr5A,Chr5B,Chr5D,Chr6A,Chr6B,Chr6D,Chr7A,Chr7B,Chr7D" \
  --output "Merged_Dataset.vcf.gz" \
  --batch-size 100 \
  --map-file "$REPO_DIR/example_sample_map.txt" \
  --non-native "OT" \
  -d &> pipeline.log

echo "Merging completed. Log written to pipeline.log"
