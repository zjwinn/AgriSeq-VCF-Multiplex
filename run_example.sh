#!/bin/bash

# ==============================================================================
# Example Wrapper Script
# This file demonstrates exactly how to call the AgriSeq_Merger.sh pipeline
# using command line flags, rather than hardcoding variables into the main script.
# ==============================================================================

# 1. Activate the required Conda environment
echo "Activating conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || eval "$(conda shell.bash hook 2>/dev/null)"
conda activate agriseq_vcf_env

# 2. Execute the pipeline using the built-in flags
# Execute the master pipeline with the provided inputs and save the log
echo "Executing AgriSeq_Merger.sh pipeline... (Output is being saved to pipeline.log)"
./AgriSeq_Merger.sh \
  --zip-dir "NAG_Zip_Files" \
  --ref-index "CHS21_OT" \
  --contigs "Chr1A,Chr1B,Chr1D,Chr2A,Chr2B,Chr2D,Chr3A,Chr3B,Chr3D,Chr4A,Chr4B,Chr4D,Chr5A,Chr5B,Chr5D,Chr6A,Chr6B,Chr6D,Chr7A,Chr7B,Chr7D" \
  --output "USDA_AgriSeq_2025_2026.vcf.gz" \
  --batch-size 100 \
  --map-file "KSM2026-GYT_Sample_Names.txt" \
  --non-native "OT" \
  -d &> pipeline.log
