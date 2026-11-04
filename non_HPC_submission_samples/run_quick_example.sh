#!/bin/bash

# ==============================================================================
# Example Wrapper Script for Fast Manifest-Driven Multiplexing
# Demonstrates calling AgriSeq_Merge_Subset_Sort_Quick.sh with command line flags.
# ==============================================================================

# Exit immediately if any command exits with a non-zero exit status
set -eo pipefail

# Determine repository root directory relative to this script location
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Verify that the primary pipeline script exists in the repository root
if [ ! -f "$REPO_DIR/AgriSeq_Merge_Subset_Sort_Quick.sh" ]; then
    
    # Print error message if pipeline script is missing
    echo "Error: AgriSeq_Merge_Subset_Sort_Quick.sh not found in $REPO_DIR."
    
    # Print advice to clone repository first
    echo "Please ensure the repository is fully cloned: git clone https://github.com/<org>/AgriSeq-VCF-Multiplex.git"
    
    # Exit with error status code 1
    exit 1

# End of pipeline script verification check
fi

# Print status message for conda environment activation
echo "Activating conda environment..."

# Source conda shell initialization profile
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || eval "$(conda shell.bash hook 2>/dev/null)"

# Activate designated conda environment
conda activate agriseq_vcf_env

# Set path to directory containing raw zip archives delivered by provider
ZIP_DIR="raw_zip_files"

# Set path to canonical 5-column target manifest (#CHROM POS ID REF ALT)
TARGETS_FILE="$REPO_DIR/target_manifests/Wheat_AgSeqv2_Targets_Manifest.txt"

# Set path to optional sample barcode to line name mapping file
MAP_FILE="$REPO_DIR/example_inputs/example_sample_map.txt"

# Set path to optional line/taxa subset list file
LINES_FILE="$REPO_DIR/example_inputs/example_lines_list.txt"

# Set path to optional headless contig sort order text file (one contig per line)
CONTIGS_FILE="$REPO_DIR/example_inputs/example_contigs_order.txt"

# Set path for final compressed and indexed output VCF
OUTPUT_VCF="Multiplexed_Quick_Biallelic.vcf.gz"

# Detect number of available CPU cores or default to 4 threads
THREADS=$(nproc 2>/dev/null || echo 4)

# Print execution notice before launching pipeline
echo "Executing AgriSeq_Merge_Subset_Sort_Quick.sh... (Output saved to quick_multiplex.log)"

# Execute fast manifest-driven multiplexing pipeline and redirect all output to log file
"$REPO_DIR/AgriSeq_Merge_Subset_Sort_Quick.sh" \
  --zip-dir "$ZIP_DIR" \
  --targets "$TARGETS_FILE" \
  --map-file "$MAP_FILE" \
  --lines "$LINES_FILE" \
  --contigs "$CONTIGS_FILE" \
  --output "$OUTPUT_VCF" \
  --threads "$THREADS" &> quick_multiplex.log

# Print success completion message
echo "Job completed successfully!"

# Print path to generated final VCF file
echo "Final VCF:   $OUTPUT_VCF"

# Print path and note for generated CSI index file
echo "CSI Index:   ${OUTPUT_VCF}.csi (Tabix .tbi omitted due to >512 Mbp wheat chromosomes)"

# Print path to companion QC audit report
echo "QC Report:   ${OUTPUT_VCF%.vcf.gz}_qc.tsv"

# Print path to generated run log file
echo "Log File:    quick_multiplex.log"
