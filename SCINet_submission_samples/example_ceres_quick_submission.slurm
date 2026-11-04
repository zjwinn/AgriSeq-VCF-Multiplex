#!/bin/bash
#SBATCH --job-name=AgriSeq_Quick
#SBATCH --output=agriseq_quick_%j.out
#SBATCH --error=agriseq_quick_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=medium
#SBATCH --account=your_project_account
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.name@usda.gov

# Print startup banner displaying hostname and current timestamp
echo "Job started on $(hostname) at $(date)"

# Check if SLURM_SUBMIT_DIR is set by SLURM scheduler
if [ -n "$SLURM_SUBMIT_DIR" ]; then
    
    # Use SLURM submit directory as repository root
    REPO_DIR="$SLURM_SUBMIT_DIR"

# Fall back to script parent directory if running interactively
else
    
    # Determine repository root relative to script location
    REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# End of repository root resolution check
fi

# Verify that the primary pipeline script exists in the repository root
if [ ! -f "$REPO_DIR/AgriSeq_Merge_Subset_Sort_Quick.sh" ]; then
    
    # Print error message indicating missing pipeline script
    echo "Error: AgriSeq_Merge_Subset_Sort_Quick.sh not found in $REPO_DIR."
    
    # Print instruction to clone repository first
    echo "Please clone the repository first: git clone https://github.com/<org>/AgriSeq-VCF-Multiplex.git"
    
    # Exit with error status code 1
    exit 1

# End of script verification check
fi

# Purge any existing loaded environment modules
module purge

# Load bcftools cluster module
module load bcftools

# Load htslib cluster module
module load htslib

# Load miniconda3 cluster module
module load miniconda3

# Source user conda initialization profile if present
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true

# Activate agriseq_vcf_env conda environment if configured
conda activate agriseq_vcf_env 2>/dev/null || true

# Set path to raw zip directory delivered by provider
ZIP_DIR="raw_zip_files"

# Set path to canonical 5-column target manifest (#CHROM POS ID REF ALT)
TARGETS_FILE="$REPO_DIR/target_manifests/Wheat_AgSeqv2_Targets_Manifest.txt"

# Set path to optional sample barcode to line name mapping file
MAP_FILE="$REPO_DIR/example_inputs/example_sample_map.txt"

# Set path to optional lines subset filter file
LINES_FILE="$REPO_DIR/example_inputs/example_lines_list.txt"

# Set path to optional headless contig sort order text file (one contig per line)
CONTIGS_FILE="$REPO_DIR/example_inputs/example_contigs_order.txt"

# Set path for final compressed and indexed output VCF
OUTPUT_VCF="Multiplexed_Quick_Biallelic.vcf.gz"

# Execute fast manifest-driven multiplexing pipeline on Ceres
"$REPO_DIR/AgriSeq_Merge_Subset_Sort_Quick.sh" \
  --zip-dir "$ZIP_DIR" \
  --targets "$TARGETS_FILE" \
  --map-file "$MAP_FILE" \
  --lines "$LINES_FILE" \
  --contigs "$CONTIGS_FILE" \
  --output "$OUTPUT_VCF" \
  --threads "$SLURM_CPUS_PER_TASK"

# Print job completion message and timestamp
echo "Job completed at $(date)"

# Print path to generated final VCF file
echo "Final VCF:   $OUTPUT_VCF"

# Print path and note for generated CSI index file
echo "CSI Index:   ${OUTPUT_VCF}.csi (Tabix .tbi omitted due to >512 Mbp wheat chromosomes)"

# Print path to companion QC audit report
echo "QC Report:   ${OUTPUT_VCF%.vcf.gz}_qc.tsv"
