#!/bin/bash

# ==============================================================================
# AgriSeq_Merge_Subset_Sort_Quick.sh
# ==============================================================================
# A fast, unified pipeline for targeted AgriSeq VCF multiplexing.
#
# PURPOSE:
#   Directly projects raw single-sample VCFs against a canonical 5-column target
#   manifest (#CHROM POS ID REF ALT) in parallel, completely bypassing the
#   multiallelic expansion, coordinate window padding, and header bloat
#   traditionally introduced by 'bcftools merge'.
#
# CORE WORKFLOW:
#   1. Scans raw single-sample VCFs (from zip archives or a directory) and
#      maps barcode IDs to line designations via an optional map file.
#   2. Subsets to requested lines if a line list is provided (or includes all).
#   3. Projects each sample against the manifest in parallel:
#      - Retains calls matching canonical REF/ALT (0/0, 0/1, 1/1).
#      - Masks off-target alleles (flanking indels, off-target SNPs) and
#        uncalled/missing loci strictly to missing data (./.).
#      - Outputs GT-only by default (ultra-fast, compact) or full FORMAT fields
#        if -f/--full-format is enabled.
#   4. Horizontally assembles the sample columns into a valid VCF 4.2 matrix.
#   5. Sorts records in user-defined or canonical contig order with monotonically
#      increasing POS coordinates.
#   6. Compresses with bgzip and indexes with CSI (explicitly omitting .tbi
#      because wheat chromosomes > 512 Mbp exceed Tabix index format limits).
#   7. Computes and exports a companion per-marker QC audit report (<output>_qc.tsv).
# ==============================================================================

# Exit immediately if any command in a pipeline exits with a non-zero status
set -eo pipefail

# Determine the absolute path to the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if the python helper script is located in the python_bin subdirectory
if [ -f "$SCRIPT_DIR/python_bin/quick_multiplex.py" ]; then
    
    # Set the path to the python helper script in python_bin
    PY_QUICK_SCRIPT="$SCRIPT_DIR/python_bin/quick_multiplex.py"

# Otherwise check if the python helper script is located in the root script directory
elif [ -f "$SCRIPT_DIR/quick_multiplex.py" ]; then
    
    # Set the path to the python helper script in the root directory
    PY_QUICK_SCRIPT="$SCRIPT_DIR/quick_multiplex.py"

# Handle case where the python helper script is not found
else
    
    # Print error message indicating missing python helper script
    echo "Error: quick_multiplex.py helper script not found in '$SCRIPT_DIR/python_bin' or '$SCRIPT_DIR'."
    
    # Exit with error status code 1
    exit 1

# End of python helper script location check
fi

# Define function to display usage instructions and options
function show_help {
    
    # Print top banner border
    echo "================================================================================"
    
    # Print banner script title
    echo " AgriSeq_Merge_Subset_Sort_Quick.sh - Fast Manifest-Driven Multiplexing"
    
    # Print banner divider border
    echo "================================================================================"
    
    # Print empty line
    echo ""
    
    # Print USAGE section header
    echo "USAGE:"
    
    # Print command call template line 1
    echo "  ./AgriSeq_Merge_Subset_Sort_Quick.sh \\"
    
    # Print required target manifest argument line
    echo "      -t <target_manifest> \\"
    
    # Print required output path argument line
    echo "      -o <output_vcf.gz> \\"
    
    # Print input archive or directory argument line
    echo "      [-z <zip_dir> | -i <input_dir>] \\"
    
    # Print optional sample mapping file argument line
    echo "      [-m <sample_map>] \\"
    
    # Print optional lines filter list argument line
    echo "      [-l <lines_list>] \\"
    
    # Print optional contigs sort order file argument line
    echo "      [-c <contigs_order_file>] \\"
    
    # Print optional full format retention argument line
    echo "      [-f] \\"
    
    # Print optional worker threads argument line
    echo "      [-T <threads>]"
    
    # Print empty line
    echo ""
    
    # Print OPTIONS section header
    echo "OPTIONS:"
    
    # Print help flag description
    echo "  -h, --help           Show this help message and exit."
    
    # Print targets flag title
    echo "  -t, --targets        (Required) Canonical 5-column target manifest"
    
    # Print targets flag format specification
    echo "                       (#CHROM POS ID REF ALT, tab-delimited),"
    
    # Print targets flag example path
    echo "                       e.g. target_manifests/Wheat_AgSeqv2_Targets_Manifest.txt."
    
    # Print output flag description
    echo "  -o, --output         (Required) Path for the final sorted output VCF (.vcf.gz)."
    
    # Print zip directory flag title
    echo "  -z, --zip-dir        Directory containing raw .zip archives delivered by provider"
    
    # Print zip directory flag requirement explanation
    echo "                       (e.g. local_work/NAG_Zip_Files). Either -z or -i is required."
    
    # Print input directory flag title
    echo "  -i, --input-dir      Directory containing pre-extracted raw single-sample .vcf.gz files."
    
    # Print input directory flag detail
    echo "                       (Can be used instead of or in addition to -z)."
    
    # Print sample map flag title
    echo "  -m, --map-file       (Optional) Two-column lookup file mapping raw VCF sample barcode IDs"
    
    # Print sample map flag structure explanation
    echo "                       to cultivar/line names (<raw_id> <line_name>)."
    
    # Print lines filter flag title
    echo "  -l, --lines          (Optional) Single-column text file of line/taxa names to subset."
    
    # Print lines filter flag behavior explanation
    echo "                       If omitted, all discovered samples are included."
    
    # Print contigs order file flag title
    echo "  -c, --contigs        (Optional) Path to a headless text file specifying contig sort order"
    
    # Print contigs order file detail line 1
    echo "                       (one contig per line). Overrides default wheat chromosome ordering."
    
    # Print full format flag title
    echo "  -f, --full-format    (Optional) Retain deep TVC FORMAT fields (GQ, DP, RO, AO, etc.)."
    
    # Print full format default mode explanation
    echo "                       Defaults to GT-only for maximum speed and TASSEL/T3 compliance."
    
    # Print worker threads flag description
    echo "  -T, --threads        (Optional) Number of parallel worker threads (default: auto-detected CPU cores)."
    
    # Print empty line
    echo ""
    
    # Print INDEXING NOTE section header
    echo "INDEXING NOTE:"
    
    # Print CSI index usage explanation line 1
    echo "  Generates a CSI index (.csi) exclusively. Tabix (.tbi) is intentionally omitted"
    
    # Print wheat chromosome size limitation explanation line 2
    echo "  because wheat chromosomes (e.g. Chr3B > 512 Mbp) exceed Tabix index coordinate limits."
    
    # Print bottom banner border
    echo "================================================================================"

# End of show_help function definition
}

# Initialize variable for target manifest file path
TARGETS_FILE=""

# Initialize variable for final output VCF file path
OUTPUT_VCF=""

# Initialize variable for raw zip archives directory
ZIP_DIR=""

# Initialize variable for pre-extracted VCF directory
INPUT_DIR=""

# Initialize variable for sample barcode to line mapping file
MAP_FILE=""

# Initialize variable for line/taxa subset filter file
LINES_FILE=""

# Initialize variable for custom contig ordering text file
CONTIGS_ORDER_FILE=""

# Initialize flag for full TVC FORMAT retention to false (GT-only default)
FULL_FORMAT=false

# Auto-detect available CPU processing cores or fall back to 4 threads
NUM_THREADS=$(nproc 2>/dev/null || echo 4)

# Loop through all command line arguments passed to the script
while [[ $# -gt 0 ]]; do
    
    # Match the argument flag
    case "$1" in
        
        # Check for help option flag
        -h|--help)
            
            # Call show_help function
            show_help
            
            # Exit script successfully with code 0
            exit 0
            
            # End help option case
            ;;
            
        # Check for target manifest option flag
        -t|--targets)
            
            # Assign manifest path from next argument
            TARGETS_FILE="$2"
            
            # Shift past flag and value
            shift 2
            
            # End targets option case
            ;;
            
        # Check for output VCF option flag
        -o|--output)
            
            # Assign output path from next argument
            OUTPUT_VCF="$2"
            
            # Shift past flag and value
            shift 2
            
            # End output option case
            ;;
            
        # Check for zip directory option flag
        -z|--zip-dir)
            
            # Assign zip directory path from next argument
            ZIP_DIR="$2"
            
            # Shift past flag and value
            shift 2
            
            # End zip directory option case
            ;;
            
        # Check for pre-extracted input directory option flag
        -i|--input-dir|--vcf-dir)
            
            # Assign input directory path from next argument
            INPUT_DIR="$2"
            
            # Shift past flag and value
            shift 2
            
            # End input directory option case
            ;;
            
        # Check for sample mapping file option flag
        -m|--map-file)
            
            # Assign map file path from next argument
            MAP_FILE="$2"
            
            # Shift past flag and value
            shift 2
            
            # End map file option case
            ;;
            
        # Check for lines subset filter option flag
        -l|--lines|--lines-file)
            
            # Assign lines file path from next argument
            LINES_FILE="$2"
            
            # Shift past flag and value
            shift 2
            
            # End lines file option case
            ;;
            
        # Check for custom contigs sort order file option flag
        -c|--contigs|--contigs-file)
            
            # Assign contigs order file path from next argument
            CONTIGS_ORDER_FILE="$2"
            
            # Shift past flag and value
            shift 2
            
            # End contigs order file option case
            ;;
            
        # Check for full format retention option flag
        -f|--full-format)
            
            # Set full format flag to true
            FULL_FORMAT=true
            
            # Shift past single flag
            shift 1
            
            # End full format option case
            ;;
            
        # Check for worker threads option flag
        -T|--threads)
            
            # Assign worker thread count from next argument
            NUM_THREADS="$2"
            
            # Shift past flag and value
            shift 2
            
            # End worker threads option case
            ;;
            
        # Handle unrecognized option arguments
        *)
            
            # Print error message with unrecognized flag
            echo "Error: Unknown parameter passed: $1"
            
            # Call show_help function
            show_help
            
            # Exit script with error code 1
            exit 1
            
            # End unrecognized option case
            ;;

    # End of case matching
    esac

# End of argument parsing loop
done

# Validate that both target manifest and output VCF arguments were provided
if [ -z "$TARGETS_FILE" ] || [ -z "$OUTPUT_VCF" ]; then
    
    # Print error message indicating missing required arguments
    echo "Error: Missing required arguments. Please specify both --targets (-t) and --output (-o)."
    
    # Call show_help function
    show_help
    
    # Exit script with error code 1
    exit 1

# End of required arguments validation check
fi

# Validate that at least one input source directory was provided
if [ -z "$ZIP_DIR" ] && [ -z "$INPUT_DIR" ]; then
    
    # Print error message indicating missing input directory
    echo "Error: Missing input directory. Please specify either --zip-dir (-z) or --input-dir (-i)."
    
    # Call show_help function
    show_help
    
    # Exit script with error code 1
    exit 1

# End of input source validation check
fi

# Validate that the specified target manifest exists on disk
if [ ! -f "$TARGETS_FILE" ]; then
    
    # Print error message indicating non-existent target manifest
    echo "Error: Target manifest '$TARGETS_FILE' does not exist."
    
    # Exit script with error code 1
    exit 1

# End of target manifest file check
fi

# Validate that the sample mapping file exists on disk if specified
if [ -n "$MAP_FILE" ] && [ ! -f "$MAP_FILE" ]; then
    
    # Print error message indicating non-existent sample map file
    echo "Error: Sample map file '$MAP_FILE' does not exist."
    
    # Exit script with error code 1
    exit 1

# End of sample map file check
fi

# Validate that the lines subset file exists on disk if specified
if [ -n "$LINES_FILE" ] && [ ! -f "$LINES_FILE" ]; then
    
    # Print error message indicating non-existent lines filter file
    echo "Error: Lines list file '$LINES_FILE' does not exist."
    
    # Exit script with error code 1
    exit 1

# End of lines filter file check
fi

# Validate that the custom contigs order file exists on disk if specified
if [ -n "$CONTIGS_ORDER_FILE" ] && [ ! -f "$CONTIGS_ORDER_FILE" ]; then
    
    # Print error message indicating non-existent contigs order file
    echo "Error: Contigs order file '$CONTIGS_ORDER_FILE' does not exist."
    
    # Exit script with error code 1
    exit 1

# End of contigs order file check
fi

# Initialize python command variable
PYTHON_CMD=""

# Check if python3 executable is available in system PATH
if command -v python3 &> /dev/null; then
    
    # Assign python3 command
    PYTHON_CMD="python3"

# Alternatively check if python executable is available in system PATH
elif command -v python &> /dev/null; then
    
    # Assign python command
    PYTHON_CMD="python"

# Handle case where no python interpreter is found
else
    
    # Print error message indicating python is missing from PATH
    echo "Error: Neither 'python3' nor 'python' was found in PATH. Please activate your conda environment."
    
    # Exit script with error code 1
    exit 1

# End of python executable check
fi

# Check if bcftools executable is available in system PATH
if ! command -v bcftools &> /dev/null; then
    
    # Print error message indicating bcftools is missing
    echo "Error: 'bcftools' was not found in PATH. Please activate your conda environment or load the module."
    
    # Exit script with error code 1
    exit 1

# End of bcftools check
fi

# Check if bgzip executable is available in system PATH
if ! command -v bgzip &> /dev/null; then
    
    # Print error message indicating bgzip is missing
    echo "Error: 'bgzip' was not found in PATH. Please activate your conda environment or load the module."
    
    # Exit script with error code 1
    exit 1

# End of bgzip check
fi

# Create a secure temporary directory managed by the operating system
TEMP_DIR=$(mktemp -d /tmp/agriseq_quick_XXXXXX)

# Register trap handler to automatically remove the temporary directory upon script exit
trap 'rm -rf "$TEMP_DIR"' EXIT

# Create parent directory for final output VCF if it does not already exist
mkdir -p "$(dirname "$OUTPUT_VCF")"

# Strip .vcf.gz suffix from output path for companion QC report file naming
QC_REPORT_FILE="${OUTPUT_VCF%.vcf.gz}"

# Strip optional .vcf suffix and append _qc.tsv extension for companion report
QC_REPORT_FILE="${QC_REPORT_FILE%.vcf}_qc.tsv"

# Print top startup banner border
echo "================================================================================"

# Print pipeline startup title
echo " Starting AgriSeq Fast Manifest-Driven Multiplexing Pipeline"

# Print current execution timestamp
echo " Date:          $(date)"

# Print target manifest file path
echo " Target Manifest: $TARGETS_FILE"

# Print final output VCF path
echo " Output VCF:      $OUTPUT_VCF"

# Print companion QC report path
echo " Companion QC:    $QC_REPORT_FILE"

# Print zip directory path or None if omitted
echo " Zip Directory:   ${ZIP_DIR:-None}"

# Print input directory path or None if omitted
echo " Input Directory: ${INPUT_DIR:-None}"

# Print sample map file path or default note
echo " Map File:        ${MAP_FILE:-None (using raw sample IDs)}"

# Print lines filter file path or default note
echo " Lines File:      ${LINES_FILE:-None (keeping all discovered lines)}"

# Print contig sort order configuration
echo " Contig Order:    ${CONTIGS_ORDER_FILE:-Default (Wheat Chr1A..Chr7D / Natural)}"

# Print active FORMAT field extraction mode
echo " FORMAT Mode:     $(if [ "$FULL_FORMAT" = true ]; then echo "Full TVC FORMAT"; else echo "GT-only"; fi)"

# Print configured worker thread count
echo " Worker Threads:  $NUM_THREADS"

# Print temporary working directory path
echo " Temp Directory:  $TEMP_DIR"

# Print bottom startup banner border
echo "================================================================================"

# Print blank line before step 1
echo ""

# Print step 1 status message
echo "[1/4] Discovering raw single-sample VCFs and indexing sample identifiers..."

# Define temporary file path for discovered samples manifest
INDEX_MANIFEST="$TEMP_DIR/discovered_samples.tsv"

# Build initial command array for index-samples subcommand
INDEX_CMD=(
    "$PYTHON_CMD" "$PY_QUICK_SCRIPT" index-samples
    --output "$INDEX_MANIFEST"
)

# Append zip-dir argument if zip directory was provided
if [ -n "$ZIP_DIR" ]; then
    
    # Append --zip-dir flag and value to command array
    INDEX_CMD+=(--zip-dir "$ZIP_DIR")

# End of zip directory argument append
fi

# Append input-dir argument if pre-extracted directory was provided
if [ -n "$INPUT_DIR" ]; then
    
    # Append --input-dir flag and value to command array
    INDEX_CMD+=(--input-dir "$INPUT_DIR")

# End of input directory argument append
fi

# Append map-file argument if sample mapping file was provided
if [ -n "$MAP_FILE" ]; then
    
    # Append --map-file flag and value to command array
    INDEX_CMD+=(--map-file "$MAP_FILE")

# End of map file argument append
fi

# Append lines-file argument if lines subset file was provided
if [ -n "$LINES_FILE" ]; then
    
    # Append --lines-file flag and value to command array
    INDEX_CMD+=(--lines-file "$LINES_FILE")

# End of lines file argument append
fi

# Execute index-samples subcommand to discover and filter samples
"${INDEX_CMD[@]}"

# Count total discovered samples in index manifest
TOTAL_SAMPLES=$(wc -l < "$INDEX_MANIFEST" | tr -d ' ')

# Validate that at least one matching sample was found
if [ "$TOTAL_SAMPLES" -eq 0 ]; then
    
    # Print error message indicating no matching samples were found
    echo "Error: Zero matching samples found in provided input archives/directories."
    
    # Exit script with error code 1
    exit 1

# End of matching sample count check
fi

# Print blank line before step 2
echo ""

# Print step 2 status message
echo "[2/4] Validating target manifest and preparing locus definitions..."

# Define temporary file path for sorted target definitions
TARGETS_CLEAN="$TEMP_DIR/targets_clean.tsv"

# Define temporary file path for ordered contig definitions
CONTIGS_FILE="$TEMP_DIR/contigs_order.txt"

# Build initial command array for prepare-manifest subcommand
MAN_CMD=(
    "$PYTHON_CMD" "$PY_QUICK_SCRIPT" prepare-manifest
    --manifest "$TARGETS_FILE"
    --output-targets "$TARGETS_CLEAN"
    --output-contigs "$CONTIGS_FILE"
)

# Append contigs-file argument if custom contig order file was provided
if [ -n "$CONTIGS_ORDER_FILE" ]; then
    
    # Append --contigs-file flag and path to command array
    MAN_CMD+=(--contigs-file "$CONTIGS_ORDER_FILE")

# End of contigs order file argument append
fi

# Execute prepare-manifest subcommand to validate and sort target loci
"${MAN_CMD[@]}"

# Count total valid targets prepared from manifest
TOTAL_TARGETS=$(wc -l < "$TARGETS_CLEAN" | tr -d ' ')

# Validate that manifest produced at least one target record
if [ "$TOTAL_TARGETS" -eq 0 ]; then
    
    # Print error message indicating empty target records
    echo "Error: Manifest contained 0 valid target records."
    
    # Exit script with error code 1
    exit 1

# End of target records check
fi

# Print blank line before step 3
echo ""

# Print step 3 status message with sample count and thread count
echo "[3/4] Projecting $TOTAL_SAMPLES samples against $TOTAL_TARGETS manifest targets in parallel ($NUM_THREADS threads)..."

# Define temporary directory path for individual sample column output files
SAMPLE_COLS_DIR="$TEMP_DIR/sample_columns"

# Create sample columns temporary directory
mkdir -p "$SAMPLE_COLS_DIR"

# Export variables required by parallel worker functions
export TARGETS_CLEAN SAMPLE_COLS_DIR FULL_FORMAT PYTHON_CMD PY_QUICK_SCRIPT

# Define worker function to project a single sample against canonical manifest targets
process_single_sample_worker() {
    
    # Assign input source type argument (ZIP or FILE)
    local source_type="$1"
    
    # Assign container archive or file path argument
    local container_path="$2"
    
    # Assign inner archive entry path argument
    local inner_path="$3"
    
    # Assign raw sample barcode ID argument
    local sample_id="$4"
    
    # Assign resolved line/cultivar name argument
    local line_name="$5"
    
    # Assign numeric sample index argument
    local sample_idx="$6"

    # Define path for single-column output file
    local col_output="$SAMPLE_COLS_DIR/col_${sample_idx}.tsv"

    # Build command array for project-sample subcommand
    local PROJ_CMD=(
        "$PYTHON_CMD" "$PY_QUICK_SCRIPT" project-sample
        --source-type "$source_type"
        --container "$container_path"
        --inner "$inner_path"
        --targets "$TARGETS_CLEAN"
        --output "$col_output"
        --line-name "$line_name"
    )

    # Check if full TVC FORMAT retention was requested
    if [ "$FULL_FORMAT" = true ]; then
        
        # Append --full-format flag to project command array
        PROJ_CMD+=(--full-format)
    
    # End of full format check
    fi

    # Execute project-sample subcommand for current sample
    "${PROJ_CMD[@]}"

# End of process_single_sample_worker function definition
}

# Export worker function for subshell access via xargs
export -f process_single_sample_worker

# Define worker function to unpack tab-delimited manifest line and invoke worker
process_manifest_line() {
    
    # Assign raw manifest line parameter
    local raw_line="$1"
    
    # Declare local variables for manifest fields
    local source_type container_path inner_path sample_id line_name sample_idx
    
    # Parse tab-delimited values into respective local variables
    IFS=$'\t' read -r source_type container_path inner_path sample_id line_name sample_idx <<< "$raw_line"
    
    # Invoke single sample worker with unpacked field parameters
    process_single_sample_worker "$source_type" "$container_path" "$inner_path" "$sample_id" "$line_name" "$sample_idx"

# End of process_manifest_line function definition
}

# Export manifest line processor function for subshell access via xargs
export -f process_manifest_line

# Attach numeric sample line index (NR) as 6th column to preserve order, ensuring inner_path is never empty
awk -F'\t' -v OFS='\t' '{inner=($3=="" ? "-" : $3); print $1, $2, inner, $4, $5, NR}' "$INDEX_MANIFEST" > "$TEMP_DIR/indexed_manifest.tsv"

# Execute parallel sample projections across worker threads using newline delimiter
xargs -d '\n' -P "$NUM_THREADS" -I {} bash -c 'process_manifest_line "$@"' _ {} < "$TEMP_DIR/indexed_manifest.tsv"

# Print blank line before step 4
echo ""

# Print step 4 status message
echo "[4/4] Assembling unified VCF matrix, applying canonical sort, and indexing..."

# Define temporary path for unsorted horizontal raw matrix
MATRIX_RAW="$TEMP_DIR/matrix_raw.tsv"

# Define temporary path for assembled unsorted VCF file
VCF_UNSORTED="$TEMP_DIR/unsorted.vcf"

# Define temporary path for canonical sorted raw VCF file
VCF_SORTED_RAW="$TEMP_DIR/sorted_raw.vcf"

# Build command array for assemble-matrix subcommand
ASM_CMD=(
    "$PYTHON_CMD" "$PY_QUICK_SCRIPT" assemble-matrix
    --targets "$TARGETS_CLEAN"
    --cols-dir "$SAMPLE_COLS_DIR"
    --total-samples "$TOTAL_SAMPLES"
    --output "$MATRIX_RAW"
)

# Check if full TVC FORMAT retention was requested
if [ "$FULL_FORMAT" = true ]; then
    
    # Append --full-format flag to assemble command array
    ASM_CMD+=(--full-format)

# End of full format check
fi

# Execute assemble-matrix subcommand to join all sample columns
"${ASM_CMD[@]}"

# Define temporary path for assembled VCF header file
VCF_HEADER="$TEMP_DIR/vcf_header.txt"

# Write base VCF 4.2 metadata and standard GT FORMAT header lines
cat << 'HEADER_EOF' > "$VCF_HEADER"
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
HEADER_EOF

# Check if full TVC FORMAT fields were requested
if [ "$FULL_FORMAT" = true ]; then

# Append full TVC FORMAT metadata header lines
cat << 'FMT_EOF' >> "$VCF_HEADER"
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=FDP,Number=1,Type=Integer,Description="Filtered Flow Depth">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference Allele Observation Count">
##FORMAT=<ID=FRO,Number=1,Type=Integer,Description="Filtered Reference Allele Observation Count">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate Allele Observation Count">
##FORMAT=<ID=FAO,Number=A,Type=Integer,Description="Filtered Alternate Allele Observation Count">
##FORMAT=<ID=SAR,Number=A,Type=Integer,Description="Alternate Allele Observation Count on Reverse Strand">
##FORMAT=<ID=SAF,Number=A,Type=Integer,Description="Alternate Allele Observation Count on Forward Strand">
##FORMAT=<ID=SRF,Number=1,Type=Integer,Description="Reference Allele Observation Count on Forward Strand">
##FORMAT=<ID=SRR,Number=1,Type=Integer,Description="Reference Allele Observation Count on Reverse Strand">
FMT_EOF

# End of full format header append
fi

# Append contig declaration header lines formatted from ordered contigs list
awk '{print "##contig=<ID=" $1 ">"}' "$CONTIGS_FILE" >> "$VCF_HEADER"

# Combine VCF header metadata with horizontal matrix data rows
cat "$VCF_HEADER" "$MATRIX_RAW" > "$VCF_UNSORTED"

# Print status message for canonical sorting
echo "  -> Sorting records into designated contig and coordinate order..."

# Execute sort-vcf subcommand passing the ordered contigs definition file
"$PYTHON_CMD" "$PY_QUICK_SCRIPT" sort-vcf \
    --input "$VCF_UNSORTED" \
    --output "$VCF_SORTED_RAW" \
    --contigs-file "$CONTIGS_FILE"

# Print status message for bgzip compression
echo "  -> Compressing with bgzip to $OUTPUT_VCF..."

# Compress sorted VCF to final destination using bgzip
bgzip -c "$VCF_SORTED_RAW" > "$OUTPUT_VCF"

# Print status message for CSI indexing
echo "  -> Creating CSI index (.csi)..."

# Generate CSI index exclusively (Tabix .tbi is omitted due to >512 Mbp wheat chromosomes)
bcftools index -f -c "$OUTPUT_VCF"

# Print status message for companion QC calculation
echo "  -> Computing companion QC audit statistics..."

# Execute generate-qc subcommand to generate companion QC TSV and console audit summary
"$PYTHON_CMD" "$PY_QUICK_SCRIPT" generate-qc \
    --vcf "$VCF_SORTED_RAW" \
    --qc-output "$QC_REPORT_FILE" \
    --final-vcf "$OUTPUT_VCF"

# Exit script successfully with code 0
exit 0
