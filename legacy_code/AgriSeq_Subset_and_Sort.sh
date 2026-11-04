#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define a function to display the help and usage section
function show_help {
    
    # Print the header line
    echo "================================================================================"
    
    # Print the script title
    echo " AgriSeq_Subset_and_Sort.sh - Non-Destructive Subsetting and Genomic Sorting"
    echo "================================================================================"
    
    # Print a blank line
    echo ""
    
    # Print the USAGE section header
    echo "USAGE:"
    
    # Print the general command structure
    echo "  ./AgriSeq_Subset_and_Sort.sh -i <input_vcf> -o <output_vcf> [-t <targets_file>] [-b] [-l <lines_file>] [-m <memory>]"
    
    # Print a blank line
    echo ""
    
    # Print the OPTIONS section header
    echo "OPTIONS:"
    
    # Print the help flag description
    echo "  -h, --help               Show this help message and exit."
    
    # Print the input VCF flag description
    echo "  -i, --input              (Required) Path to input VCF file (.vcf or .vcf.gz)."
    
    # Print the output VCF flag description
    echo "  -o, --output             (Required) Path for final output VCF file (.vcf.gz or .vcf)."
    
    # Print the targets list flag description
    echo "  -t, --targets            (Optional) Tab-delimited target manifest dynamically inspected for"
    echo "                           required columns (#CHROM POS ID by default; #CHROM POS ID REF ALT"
    echo "                           when -b/--biallelic-targets is enabled). Single-column ID lists"
    echo "                           are disallowed due to off-target amplicon indels."
    echo "                           (e.g. target_manifests/Wheat_AgSeqv2_Targets_Manifest.txt or"
    echo "                           example_inputs/example_target_manifest.txt)."
    
    # Print the biallelic standardization flag description
    echo "  -b, --biallelic-targets  (Optional) Enforce strict canonical biallelic REF and ALT alleles"
    echo "                           defined in the target manifest (--targets required). Trims merge-padded"
    echo "                           REF/ALT windows, strips off-target secondary alleles, slices Number=A/R/G"
    echo "                           FORMAT/INFO fields to the target ALT, masks sample genotypes carrying"
    echo "                           off-target alleles completely to ./., and writes a companion QC TSV"
    echo "                           (<output_prefix>_biallelic_qc.tsv) for T3/Breedbase, GWAS, and GS."
    
    # Print the lines list flag description
    echo "  -l, --lines              (Optional) Single-column text file of line/taxa names"
    echo "                           (e.g. example_inputs/example_lines_list.txt)."
    
    # Print the memory flag description
    echo "  -m, --memory             (Optional) Memory limit for bcftools sorting buffer (default: 16g)."
    
    # Print a blank line
    echo ""
    
    # Print the PYTHON SCRIPT USAGE header
    echo "PYTHON SCRIPT USAGE (standardize_biallelic_targets.py):"
    
    # Print explanation for internal and standalone python helper execution
    echo "  When -t/--targets is provided, the pipeline calls python_bin/standardize_biallelic_targets.py"
    echo "  to verify coordinate + marker ID pairing (and enforce biallelic REF/ALT when -b is passed)."
    echo "  View standalone options with: python python_bin/standardize_biallelic_targets.py -h"
    
    # Print a blank line
    echo ""
    
    # Print the DETAILS header
    echo "PIPELINE DETAILS & INTEGRITY CHECKS:"
    
    # Print description of dynamic manifest validation
    echo "  1. Dynamic Manifest & Coordinate+ID Validation: Dynamically inspects TSV column headers"
    echo "     (#CHROM, POS, ID, and REF/ALT when -b is active), verifies positive integer coordinates,"
    echo "     validates contigs against the VCF header, and enforces both coordinate and marker ID"
    echo "     pairing to reject adjacent off-target amplicon indels."
    
    # Print description of sample filtering
    echo "  2. Sample Line Validation: Validates single-column line files, collapses duplicate"
    echo "     taxa, and intersects sample names with the VCF header, reporting missing samples."
    
    # Print description of TASSEL compatibility
    echo "  3. Genomic Ordering: Sorts records in strict physical coordinate order with"
    echo "     natural contig sequencing so VCFs open cleanly in the TASSEL GUI."
    
    # Print description of preservation vs biallelic standardization
    echo "  4. Non-Destructive Default vs. Strict Biallelic Mode (-b): Preserves 100% of multiallelic"
    echo "     calls and headers by default, or standardizes loci strictly to canonical target REF/ALT"
    echo "     alleles (masking off-target calls to ./.) when -b/--biallelic-targets is specified."
    
    # Print a blank line
    echo ""
    
    # Print the REQUIREMENTS header
    echo "REQUIREMENTS:"
    
    # Print the environment requirement description
    echo "  Must be run within an environment containing 'bcftools', 'htslib', and 'python' (via Conda or module load)."
    
    # Print the final border line
    echo "================================================================================"
    
# End of the show_help function
}

# Initialize the INPUT_VCF variable
INPUT_VCF=""

# Initialize the OUTPUT_VCF variable
OUTPUT_VCF=""

# Initialize the TARGETS_FILE variable
TARGETS_FILE=""

# Initialize the LINES_FILE variable
LINES_FILE=""

# Initialize the SORT_MEMORY variable with a default of 16g
SORT_MEMORY="16g"

# Initialize the BIALLELIC_TARGETS flag with a default of false (non-destructive mode)
BIALLELIC_TARGETS=false

# Initialize the companion QC report path variable
QC_REPORT_FILE=""

# Loop through all command line arguments
while [[ "$#" -gt 0 ]]; do
    
    # Check the current flag
    case $1 in
        
        # If the help flag is passed
        -h|--help)
            
            # Show the help section
            show_help
            
            # Exit the script successfully
            exit 0
            
            # End of help case
            ;;
            
        # If the input VCF flag is passed
        -i|--input)
            
            # Assign the next argument to INPUT_VCF
            INPUT_VCF="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of input VCF case
            ;;
            
        # If the output VCF flag is passed
        -o|--output)
            
            # Assign the next argument to OUTPUT_VCF
            OUTPUT_VCF="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of output VCF case
            ;;
            
        # If the targets list flag is passed
        -t|--targets)
            
            # Assign the next argument to TARGETS_FILE
            TARGETS_FILE="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of targets list case
            ;;
            
        # If the strict biallelic standardization flag is passed
        -b|--biallelic-targets)
            
            # Enable strict biallelic target standardization
            BIALLELIC_TARGETS=true
            
            # Shift the arguments by 1
            shift 1
            
            # End of biallelic targets case
            ;;
            
        # If the lines list flag is passed
        -l|--lines)
            
            # Assign the next argument to LINES_FILE
            LINES_FILE="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of lines list case
            ;;
            
        # If the memory flag is passed
        -m|--memory)
            
            # Assign the next argument to SORT_MEMORY
            SORT_MEMORY="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of memory case
            ;;
            
        # If an unknown flag is passed
        *)
            
            # Print an error message about the unknown parameter
            echo "Error: Unknown parameter passed: $1"
            
            # Show the help section
            show_help
            
            # Exit the script with an error code
            exit 1
            
            # End of unknown flag case
            ;;
            
    # End of the case statement
    esac
    
# End of the argument parsing while loop
done

# ==============================================================================
# PRE-FLIGHT SANITY CHECKS
# ==============================================================================

# Sanity Check 1: Check if required arguments are provided
if [ -z "$INPUT_VCF" ] || [ -z "$OUTPUT_VCF" ]; then
    
    # Print an error indicating missing arguments
    echo "Error: Missing required arguments (--input and --output)."
    
    # Show the help section to guide the user
    show_help
    
    # Exit the script with an error code
    exit 1
    
# End of required arguments check
fi

# Sanity Check 2: Ensure the input VCF file actually exists
if [ ! -f "$INPUT_VCF" ]; then
    
    # Print an error message
    echo "Error: Input VCF file '$INPUT_VCF' does not exist."
    
    # Exit the script with an error code
    exit 1
    
# End of input file existence check
fi

# Sanity Check 3: Ensure the input VCF file is not empty
if [ ! -s "$INPUT_VCF" ]; then
    
    # Print an error message
    echo "Error: Input VCF file '$INPUT_VCF' is empty (0 bytes)."
    
    # Exit the script with an error code
    exit 1
    
# End of input file non-empty check
fi

# Sanity Check 4: Ensure input and output VCF paths do not point to the exact same file
if [ -e "$OUTPUT_VCF" ] && [ "$INPUT_VCF" -ef "$OUTPUT_VCF" ]; then
    
    # Print an error message
    echo "Error: Input VCF and Output VCF point to the exact same file path."
    
    # Explain why this is not permitted
    echo "Cannot overwrite the input VCF file in place while streaming records."
    
    # Exit the script with an error code
    exit 1
    
# End of input/output path collision check
fi

# Sanity Check 5: Ensure the output directory exists or can be created with write permission
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")

# Check if parent directory exists; if not, create it
if [ ! -d "$OUTPUT_DIR" ]; then
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR" 2>/dev/null || {
        
        # Print error if directory creation fails
        echo "Error: Unable to create output directory '$OUTPUT_DIR'. Check filesystem permissions."
        
        # Exit with error code
        exit 1
    }
    
# End of directory creation check
fi

# Check if the output directory is writable
if [ ! -w "$OUTPUT_DIR" ]; then
    
    # Print error message
    echo "Error: Output directory '$OUTPUT_DIR' is not writable."
    
    # Exit the script with an error code
    exit 1
    
# End of output directory write permission check
fi

# Sanity Check 6: Ensure --biallelic-targets (-b) is accompanied by --targets (-t)
if [ "$BIALLELIC_TARGETS" = true ] && [ -z "$TARGETS_FILE" ]; then
    
    # Print an error message explaining the dependency
    echo "Error: --biallelic-targets (-b) requires a target manifest (--targets / -t) containing canonical REF and ALT columns."
    
    # Exit the script with an error code
    exit 1
    
# End of biallelic flag dependency check
fi

# Sanity Check 7: If a targets file is provided, ensure it exists, is non-empty, and Python helper is available
if [ -n "$TARGETS_FILE" ]; then
    
    # Check if target file exists
    if [ ! -f "$TARGETS_FILE" ]; then
        
        # Print an error message
        echo "Error: Targets file '$TARGETS_FILE' does not exist."
        
        # Exit the script with an error code
        exit 1
        
    # End of target file existence check
    fi
    
    # Check if target file is empty
    if [ ! -s "$TARGETS_FILE" ]; then
        
        # Print an error message
        echo "Error: Targets file '$TARGETS_FILE' is empty (0 bytes)."
        
        # Exit the script with an error code
        exit 1
        
    # End of target file empty check
    fi
    
    # Resolve Python interpreter for coordinate + ID verification and biallelic standardization
    if command -v python3 &> /dev/null; then
        PYTHON_CMD="python3"
    elif command -v python &> /dev/null; then
        PYTHON_CMD="python"
    else
        echo "Error: 'python' (or 'python3') is not available. Please activate your conda environment or load python."
        exit 1
    fi
    
    # Resolve directory containing this script
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    
    # Locate standardize_biallelic_targets.py helper script (check python_bin/ first, then script root)
    if [ -f "$SCRIPT_DIR/python_bin/standardize_biallelic_targets.py" ]; then
        PY_BIALLELIC_SCRIPT="$SCRIPT_DIR/python_bin/standardize_biallelic_targets.py"
    elif [ -f "$SCRIPT_DIR/standardize_biallelic_targets.py" ]; then
        PY_BIALLELIC_SCRIPT="$SCRIPT_DIR/standardize_biallelic_targets.py"
    else
        echo "Error: standardize_biallelic_targets.py not found in '$SCRIPT_DIR/python_bin' or '$SCRIPT_DIR'."
        exit 1
    fi
    
# End of targets file parameter check
fi

# Sanity Check 7: If a lines file is provided, ensure it exists and is non-empty
if [ -n "$LINES_FILE" ]; then
    
    # Check if lines file exists
    if [ ! -f "$LINES_FILE" ]; then
        
        # Print an error message
        echo "Error: Lines file '$LINES_FILE' does not exist."
        
        # Exit the script with an error code
        exit 1
        
    # End of lines file existence check
    fi
    
    # Check if lines file is empty
    if [ ! -s "$LINES_FILE" ]; then
        
        # Print an error message
        echo "Error: Lines file '$LINES_FILE' is empty (0 bytes)."
        
        # Exit the script with an error code
        exit 1
        
    # End of lines file empty check
    fi
    
# End of lines file parameter check
fi

# Sanity Check 8: Ensure bcftools is installed and available in PATH
if ! command -v bcftools &> /dev/null; then
    
    # Print an error message
    echo "Error: 'bcftools' is not available. Please activate your conda environment or load the module (module load bcftools)."
    
    # Exit the script with an error code
    exit 1
    
# End of bcftools check
fi

# Sanity Check 9: Ensure bgzip is installed and available in PATH
if ! command -v bgzip &> /dev/null; then
    
    # Print an error message
    echo "Error: 'bgzip' is not available. Please activate your conda environment or load the module (module load htslib)."
    
    # Exit the script with an error code
    exit 1
    
# End of bgzip check
fi

# Sanity Check 10: Ensure the input file is a valid, readable VCF
if ! bcftools view -h "$INPUT_VCF" > /dev/null 2>&1; then
    
    # Print an error message
    echo "Error: Input file '$INPUT_VCF' is not a valid VCF file or the header is corrupted."
    
    # Print guidance hint
    echo "Please verify that the file is an intact VCF (.vcf) or block-gzipped VCF (.vcf.gz)."
    
    # Exit the script with an error code
    exit 1
    
# End of VCF header validity check
fi

# ==============================================================================
# PIPELINE EXECUTION
# ==============================================================================

# Create a secure temporary directory managed by the OS
TEMP_DIR=$(mktemp -d)

# Set a trap to ensure the temporary directory is safely deleted whenever the script exits
trap 'rm -rf "$TEMP_DIR"' EXIT

# Print the startup banner
echo "================================================================================"
echo " Starting AgriSeq Subsetting & Sorting Pipeline"
echo " Date:          $(date)"
echo " Input VCF:     $INPUT_VCF"
echo " Output VCF:    $OUTPUT_VCF"
echo " Targets File:  ${TARGETS_FILE:-None (keeping all loci)}"
echo " Biallelic Mode:$BIALLELIC_TARGETS"
echo " Lines File:    ${LINES_FILE:-None (keeping all sample lines)}"
echo " Sort Memory:   $SORT_MEMORY"
echo " Temp Directory:$TEMP_DIR"
echo "================================================================================"

# Print section 1 header
echo "[1/4] Preparing filter parameters and validating manifest and line files..."

# Initialize array for bcftools filter options
FILTER_OPTS=()

# Check if a targets file was provided
if [ -n "$TARGETS_FILE" ]; then
    
    # Clean targets file into a temporary list stripping Windows carriage returns, UTF-8 BOM, and empty lines
    sed '1s/^\xEF\xBB\xBF//' "$TARGETS_FILE" | tr -d '\r' | grep -ve '^\s*$' > "$TEMP_DIR/cleaned_targets.txt" || true
    
    # Strip leading ## metadata comment lines if present
    awk '!/^##/' "$TEMP_DIR/cleaned_targets.txt" > "$TEMP_DIR/manifest_body.txt"
    
    # Ensure at least one non-metadata line exists in the manifest
    if [ ! -s "$TEMP_DIR/manifest_body.txt" ]; then
        
        # Print error message
        echo "Error: Target manifest '$TARGETS_FILE' contains no header or target records."
        
        # Exit with error code
        exit 1
        
    # End of non-empty body check
    fi
    
    # Inspect the first non-metadata row to dynamically detect delimiter and column headers
    FIRST_MANIFEST_LINE=$(head -n 1 "$TEMP_DIR/manifest_body.txt")
    
    # Use strict tab delimiter when tabs are present; otherwise allow general whitespace
    AWK_FS_ARGS=()
    if [[ "$FIRST_MANIFEST_LINE" == *$'\t'* ]]; then
        AWK_FS_ARGS=(-F'\t')
    fi
    
    # Dynamically inspect header presence and resolve 1-based column indices (CHROM, POS, ID, REF, ALT)
    HEADER_PARSE=$(awk "${AWK_FS_ARGS[@]}" 'NR==1 {
        has_hdr = ($0 ~ /^#/) ? 1 : 0;
        chrom_c = 0; pos_c = 0; id_c = 0; ref_c = 0; alt_c = 0;
        for (i = 1; i <= NF; i++) {
            col = $i;
            gsub(/^[ \t#]+|[ \t]+$/, "", col);
            u = toupper(col);
            if ((u == "CHROM" || u == "CHR") && chrom_c == 0) chrom_c = i;
            else if (u == "POS" && pos_c == 0) pos_c = i;
            else if (u == "ID" && id_c == 0) id_c = i;
            else if (u == "REF" && ref_c == 0) ref_c = i;
            else if (u == "ALT" && alt_c == 0) alt_c = i;
        }
        if (chrom_c > 0 && pos_c > 0) has_hdr = 1;
        print has_hdr, NF, chrom_c, pos_c, id_c, ref_c, alt_c;
        exit;
    }' "$TEMP_DIR/manifest_body.txt")
    
    # Unpack parsed header metadata variables
    read -r HAS_MANIFEST_HEADER FIRST_NF CHROM_COL POS_COL ID_COL REF_COL ALT_COL <<< "$HEADER_PARSE"
    
    # Reject single-column or 2-column marker lists with an informative explanation
    if [ "$FIRST_NF" -lt 3 ]; then
        
        # Print descriptive error banner explaining why coordinate + ID columns are mandatory
        echo "================================================================================"
        echo "Error: Target manifest format violation in '$TARGETS_FILE'!"
        echo "Found $FIRST_NF column(s) on first line: '$FIRST_MANIFEST_LINE'"
        echo ""
        echo "The target manifest MUST contain at least 3 tab-delimited columns:"
        echo "  Default mode   : #CHROM  POS  ID  [REF  ALT ...]"
        echo "  Biallelic (-b) : #CHROM  POS  ID  REF  ALT"
        echo ""
        echo "WHY COORDINATES + ID ARE REQUIRED:"
        echo "Passing a single-column marker name list is strictly prohibited because"
        echo "Thermo Fisher TVC VCFs frequently report multiple off-target amplicon indels"
        echo "at the same locus that share the exact same marker ID as the target SNP."
        echo "Filtering by both physical coordinates (#CHROM, POS) and ID is mandatory to"
        echo "isolate the true target variant and reject adjacent amplicon indels."
        echo ""
        echo "Example valid manifest lines:"
        echo "  #CHROM             POS      ID                 REF  ALT"
        echo "  Chr1A              1169383  CHS21_1A001169383  C    T"
        echo "  TEL10_7E739933655  255      TEL10_7E739933655  C    T"
        echo "================================================================================"
        
        # Exit with error code
        exit 1
        
    # End of <3 column check
    fi
    
    # Validate required columns based on header presence and biallelic mode (-b)
    if [ "$HAS_MANIFEST_HEADER" -eq 1 ]; then
        
        # Verify #CHROM/CHROM, POS, and ID columns were found in the header
        if [ "$CHROM_COL" -eq 0 ] || [ "$POS_COL" -eq 0 ] || [ "$ID_COL" -eq 0 ]; then
            echo "================================================================================"
            echo "Error: Target manifest '$TARGETS_FILE' is missing required header column(s)!"
            echo "Detected header: '$FIRST_MANIFEST_LINE'"
            echo "Required columns: #CHROM (or CHROM), POS, and ID (plus REF and ALT when using -b)."
            echo "================================================================================"
            exit 1
        fi
        
        # Verify REF and ALT columns exist when -b/--biallelic-targets is enabled
        if [ "$BIALLELIC_TARGETS" = true ] && { [ "$REF_COL" -eq 0 ] || [ "$ALT_COL" -eq 0 ]; }; then
            echo "================================================================================"
            echo "Error: --biallelic-targets (-b) requires canonical REF and ALT columns in '$TARGETS_FILE'!"
            echo "Detected header: '$FIRST_MANIFEST_LINE'"
            echo ""
            echo "Please provide a 5-column canonical target manifest (#CHROM POS ID REF ALT),"
            echo "such as target_manifests/Wheat_AgSeqv2_Targets_Manifest.txt, or omit -b to run"
            echo "default non-destructive coordinate + ID subsetting."
            echo "================================================================================"
            exit 1
        fi
        
        # Extract data rows following the header row, skipping any inline # comment lines
        tail -n +2 "$TEMP_DIR/manifest_body.txt" | awk '!/^#/' > "$TEMP_DIR/manifest_data.txt"
        
    else
        
        # Headerless manifest fallback: enforce >= 5 columns when -b is active
        if [ "$BIALLELIC_TARGETS" = true ] && [ "$FIRST_NF" -lt 5 ]; then
            echo "================================================================================"
            echo "Error: --biallelic-targets (-b) requires 5 columns (CHROM POS ID REF ALT) in"
            echo "headerless target manifest '$TARGETS_FILE', but found $FIRST_NF column(s)!"
            echo "================================================================================"
            exit 1
        fi
        
        # Assign positional column indices for headerless manifest
        CHROM_COL=1
        POS_COL=2
        ID_COL=3
        REF_COL=4
        ALT_COL=5
        
        # Extract non-comment data rows
        awk '!/^#/' "$TEMP_DIR/manifest_body.txt" > "$TEMP_DIR/manifest_data.txt"
        
    # End of header vs headerless column resolution
    fi
    
    # Count total non-comment target data rows
    TARGET_COUNT=$(wc -l < "$TEMP_DIR/manifest_data.txt" | tr -d ' ')
    
    # Sanity Check: Ensure target manifest has at least one valid data row
    if [ "$TARGET_COUNT" -eq 0 ]; then
        
        # Print error message
        echo "Error: Target manifest '$TARGETS_FILE' contains no target data rows."
        
        # Print guidance hint
        echo "Please provide a valid target manifest containing '#CHROM POS ID [REF ALT]' rows."
        
        # Exit with error code
        exit 1
        
    # End of target count check
    fi
    
    # Determine minimum required column index across every data row
    MIN_REQ_COLS="$CHROM_COL"
    for c_idx in "$POS_COL" "$ID_COL"; do
        if [ "$c_idx" -gt "$MIN_REQ_COLS" ]; then MIN_REQ_COLS="$c_idx"; fi
    done
    if [ "$BIALLELIC_TARGETS" = true ]; then
        for c_idx in "$REF_COL" "$ALT_COL"; do
            if [ "$c_idx" -gt "$MIN_REQ_COLS" ]; then MIN_REQ_COLS="$c_idx"; fi
        done
    fi
    
    # Sanity Check: Verify every data row contains at least MIN_REQ_COLS columns
    INVALID_COL_LINE=$(awk "${AWK_FS_ARGS[@]}" -v min_nf="$MIN_REQ_COLS" '{if (NF < min_nf) {print NR, NF, $0; exit}}' "$TEMP_DIR/manifest_data.txt")
    
    # If any data row has fewer columns than required
    if [ -n "$INVALID_COL_LINE" ]; then
        BAD_LINE_NUM=$(echo "$INVALID_COL_LINE" | awk '{print $1}')
        FOUND_COLS=$(echo "$INVALID_COL_LINE" | awk '{print $2}')
        BAD_LINE_TEXT=$(echo "$INVALID_COL_LINE" | cut -d' ' -f3-)
        echo "Error: Target manifest row $BAD_LINE_NUM in '$TARGETS_FILE' has $FOUND_COLS column(s), but at least $MIN_REQ_COLS are required: '$BAD_LINE_TEXT'"
        exit 1
    fi
    
    # Sanity Check: POS column must be a positive integer on every data row
    NON_NUMERIC_POS=$(awk "${AWK_FS_ARGS[@]}" -v pc="$POS_COL" '{val=$pc; gsub(/^[ \t]+|[ \t]+$/, "", val); if (val !~ /^[0-9]+$/ || val <= 0) {print NR, val; exit}}' "$TEMP_DIR/manifest_data.txt")
    
    # If a non-numeric or non-positive position was detected
    if [ -n "$NON_NUMERIC_POS" ]; then
        
        # Extract the data row number and invalid position value
        BAD_POS_LINE=$(echo "$NON_NUMERIC_POS" | awk '{print $1}')
        BAD_POS_VAL=$(echo "$NON_NUMERIC_POS" | awk '{print $2}')
        
        # Print error message
        echo "Error: Invalid physical position '$BAD_POS_VAL' at data row $BAD_POS_LINE in '$TARGETS_FILE'."
        echo "Column $POS_COL (POS) must be a positive integer genomic coordinate."
        
        # Exit with error code
        exit 1
        
    # End of POS numeric check
    fi
    
    # Sanity Check: Detect duplicate (CHROM, POS) coordinates in the manifest
    DUP_COORDS=$(awk "${AWK_FS_ARGS[@]}" -v cc="$CHROM_COL" -v pc="$POS_COL" '{c=$cc; p=$pc; gsub(/^[ \t]+|[ \t]+$/, "", c); gsub(/^[ \t]+|[ \t]+$/, "", p); print c, p}' "$TEMP_DIR/manifest_data.txt" | LC_ALL=C sort | uniq -c | awk '$1 > 1 {print $2 ":" $3 " (" $1 " entries)"}')
    
    # If duplicate coordinates exist in the manifest
    if [ -n "$DUP_COORDS" ]; then
        
        # Count the number of duplicate coordinate pairs
        DUP_COUNT=$(echo "$DUP_COORDS" | wc -l | tr -d ' ')
        
        # Print a warning notice
        echo "  -> Notice: Target manifest contains $DUP_COUNT duplicate coordinate pairs (CHROM:POS):"
        
        # Print first few duplicate coordinates
        echo "$DUP_COORDS" | head -n 5 | while read -r d; do echo "       $d"; done
        
        # If there are more than 5 duplicates, note the remainder
        if [ "$DUP_COUNT" -gt 5 ]; then
            echo "       ... and $(($DUP_COUNT - 5)) additional duplicate coordinates."
        fi
        
    # End of duplicate coordinates check
    fi
    
    # Sanity Check: Verify target manifest chromosomes exist in the input VCF header
    awk "${AWK_FS_ARGS[@]}" -v cc="$CHROM_COL" '{c=$cc; gsub(/^[ \t]+|[ \t]+$/, "", c); print c}' "$TEMP_DIR/manifest_data.txt" | LC_ALL=C sort -u > "$TEMP_DIR/manifest_chroms.txt"
    TOTAL_MANIFEST_CHROMS=$(wc -l < "$TEMP_DIR/manifest_chroms.txt" | tr -d ' ')
    
    # Extract contig IDs from input VCF header
    bcftools view -h "$INPUT_VCF" | grep '^##contig=' | sed -E 's/.*ID=([^,>]+).*/\1/' | LC_ALL=C sort -u > "$TEMP_DIR/vcf_contigs.txt" || true
    TOTAL_VCF_CONTIGS=$(wc -l < "$TEMP_DIR/vcf_contigs.txt" | tr -d ' ')
    
    # If contigs are declared in the VCF header, check for overlap
    if [ "$TOTAL_VCF_CONTIGS" -gt 0 ]; then
        
        # Find contigs that match between manifest and VCF header
        awk 'FNR==NR {vcf[$1]=1; next} {if ($1 in vcf) print $1}' "$TEMP_DIR/vcf_contigs.txt" "$TEMP_DIR/manifest_chroms.txt" > "$TEMP_DIR/matched_chroms.txt"
        MATCHED_CHROM_COUNT=$(wc -l < "$TEMP_DIR/matched_chroms.txt" | tr -d ' ')
        
        # If zero chromosomes match, throw fatal error
        if [ "$MATCHED_CHROM_COUNT" -eq 0 ]; then
            
            # Print coordinate mismatch error banner
            echo "================================================================================"
            echo "Error: Reference chromosome coordinate mismatch!"
            echo "None of the $TOTAL_MANIFEST_CHROMS chromosome/contig names in '$TARGETS_FILE' match"
            echo "the contigs declared in the header of '$INPUT_VCF'."
            echo ""
            echo "Example contig IDs in input VCF header:"
            head -n 5 "$TEMP_DIR/vcf_contigs.txt" | while read -r c; do echo "  - $c"; done
            echo ""
            echo "Example chromosome names in target manifest:"
            head -n 5 "$TEMP_DIR/manifest_chroms.txt" | while read -r c; do echo "  - $c"; done
            echo ""
            echo "Please ensure manifest contig names match the reference FASTA contig IDs."
            echo "================================================================================"
            
            # Exit with error code
            exit 1
            
        # If only a subset of manifest chromosomes match VCF contigs
        elif [ "$MATCHED_CHROM_COUNT" -lt "$TOTAL_MANIFEST_CHROMS" ]; then
            
            # Calculate number of unmatched contigs
            UNMATCHED_COUNT=$((TOTAL_MANIFEST_CHROMS - MATCHED_CHROM_COUNT))
            
            # Print informational notice
            echo "  -> Notice: $UNMATCHED_COUNT of $TOTAL_MANIFEST_CHROMS manifest contigs are not declared in VCF header."
            
        # End of chromosome match count check
        fi
        
    # End of VCF contig presence check
    fi
    
    # Print status message confirming dynamically validated target manifest
    if [ "$BIALLELIC_TARGETS" = true ]; then
        echo "  -> Validated $FIRST_NF-column target manifest (strict biallelic REF/ALT mode): $TARGET_COUNT loci across $TOTAL_MANIFEST_CHROMS contigs."
    else
        echo "  -> Validated $FIRST_NF-column target manifest (non-destructive coordinate + ID mode): $TARGET_COUNT loci across $TOTAL_MANIFEST_CHROMS contigs."
    fi
    
    # Build 1-based inclusive [POS-1, POS+1] region windows for fast C-level bcftools pre-filtering
    # (--targets-overlap 1 retains exact coordinates, 5'-anchored blocks spanning POS, and 1-bp repeat shifts
    #  before standardize_biallelic_targets.py resolves exact coordinate + ID matches)
    awk "${AWK_FS_ARGS[@]}" -v OFS='\t' -v cc="$CHROM_COL" -v pc="$POS_COL" '{
        c = $cc; p = $pc;
        gsub(/^[ \t]+|[ \t]+$/, "", c);
        gsub(/^[ \t]+|[ \t]+$/, "", p);
        p_num = p + 0;
        start_p = (p_num > 1) ? (p_num - 1) : 1;
        end_p = p_num + 1;
        print c, start_p, end_p;
    }' "$TEMP_DIR/manifest_data.txt" > "$TEMP_DIR/target_regions.tsv"
    
    # Add target region pre-filter flags (-T with overlap mode 1)
    FILTER_OPTS+=(-T "$TEMP_DIR/target_regions.tsv" --targets-overlap 1)
    
# End of targets file processing
fi

# Check if a lines file was provided
if [ -n "$LINES_FILE" ]; then
    
    # Clean line list into a temporary file stripping Windows carriage returns and empty lines
    grep -ve '^\s*$' "$LINES_FILE" | tr -d '\r' > "$TEMP_DIR/cleaned_lines.txt"
    
    # Count the number of lines requested
    REQUESTED_LINES=$(wc -l < "$TEMP_DIR/cleaned_lines.txt" | tr -d ' ')
    
    # Sanity Check: Ensure lines file has at least one valid line entry
    if [ "$REQUESTED_LINES" -eq 0 ]; then
        
        # Print error message
        echo "Error: Lines file '$LINES_FILE' contains no sample entries (only whitespace or empty lines)."
        
        # Exit with error code
        exit 1
        
    # End of requested lines count check
    fi
    
    # Sanity Check: If lines file has multiple whitespace columns, isolate column 1
    MULTI_COL_LINES=$(awk '{if (NF > 1) {print NR, NF; exit}}' "$TEMP_DIR/cleaned_lines.txt")
    
    # If multiple columns were found on a line
    if [ -n "$MULTI_COL_LINES" ]; then
        
        # Print notice to user
        echo "  -> Notice: Lines file '$LINES_FILE' contains multiple columns; extracting column 1 as sample IDs."
        
        # Extract column 1
        awk '{print $1}' "$TEMP_DIR/cleaned_lines.txt" > "$TEMP_DIR/cleaned_lines_col1.txt"
        mv "$TEMP_DIR/cleaned_lines_col1.txt" "$TEMP_DIR/cleaned_lines.txt"
        
    # End of multi-column lines check
    fi
    
    # Sanity Check: Detect and collapse duplicate sample names in lines file
    DUP_SAMPLES=$(LC_ALL=C sort "$TEMP_DIR/cleaned_lines.txt" | uniq -d)
    
    # If duplicate sample designations exist
    if [ -n "$DUP_SAMPLES" ]; then
        
        # Count the number of duplicate lines
        DUP_SAMPLE_COUNT=$(echo "$DUP_SAMPLES" | wc -l | tr -d ' ')
        
        # Print notice to user
        echo "  -> Notice: Collapsing $DUP_SAMPLE_COUNT duplicate sample designations in lines file."
        
        # Collapse duplicates into unique list
        LC_ALL=C sort -u "$TEMP_DIR/cleaned_lines.txt" > "$TEMP_DIR/unique_lines.txt"
        mv "$TEMP_DIR/unique_lines.txt" "$TEMP_DIR/cleaned_lines.txt"
        
        # Recalculate requested lines count
        REQUESTED_LINES=$(wc -l < "$TEMP_DIR/cleaned_lines.txt" | tr -d ' ')
        
    # End of duplicate sample names check
    fi
    
    # Extract all sample names currently present in the input VCF
    bcftools query -l "$INPUT_VCF" > "$TEMP_DIR/all_vcf_samples.txt"
    TOTAL_VCF_SAMPLES=$(wc -l < "$TEMP_DIR/all_vcf_samples.txt" | tr -d ' ')
    
    # Intersect requested lines with samples actually present in the input VCF
    awk 'FNR==NR { vcf[$1]=1; next } { if ($1 in vcf) print $1 }' "$TEMP_DIR/all_vcf_samples.txt" "$TEMP_DIR/cleaned_lines.txt" > "$TEMP_DIR/matched_samples.txt"
    
    # Count how many requested lines matched
    MATCHED_COUNT=$(wc -l < "$TEMP_DIR/matched_samples.txt" | tr -d ' ')
    
    # Sanity Check: If zero lines matched, throw fatal error
    if [ "$MATCHED_COUNT" -eq 0 ]; then
        
        # Print error banner
        echo "================================================================================"
        echo "Error: None of the $REQUESTED_LINES lines in '$LINES_FILE' matched samples in '$INPUT_VCF'."
        echo ""
        echo "Example sample names present in VCF header ($TOTAL_VCF_SAMPLES total):"
        head -n 5 "$TEMP_DIR/all_vcf_samples.txt" | while read -r s; do echo "  - $s"; done
        echo ""
        echo "Example sample names requested in lines file:"
        head -n 5 "$TEMP_DIR/cleaned_lines.txt" | while read -r s; do echo "  - $s"; done
        echo ""
        echo "Please verify that the line designations match the sample names in the VCF header."
        echo "================================================================================"
        
        # Exit with error code
        exit 1
        
    # End of zero matched count check
    fi
    
    # If some requested lines were missing from the VCF, print a warning notice
    if [ "$MATCHED_COUNT" -lt "$REQUESTED_LINES" ]; then
        
        # Calculate number of missing lines
        MISSING_COUNT=$((REQUESTED_LINES - MATCHED_COUNT))
        
        # Print warning notice
        echo "  -> Warning: $MISSING_COUNT of $REQUESTED_LINES requested lines were NOT found in input VCF."
        echo "     Example missing sample lines:"
        awk 'FNR==NR { vcf[$1]=1; next } { if (!($1 in vcf)) print "       - "$1 }' "$TEMP_DIR/all_vcf_samples.txt" "$TEMP_DIR/cleaned_lines.txt" | head -n 5
        
    # End of partial match check
    fi
    
    # Print status message
    echo "  -> Subsetting $MATCHED_COUNT of $REQUESTED_LINES requested lines present in VCF..."
    
    # Add sample filtering flag to bcftools arguments
    FILTER_OPTS+=(-S "$TEMP_DIR/matched_samples.txt" --force-samples)
    
# End of lines file processing
fi

# Print section 2 header
echo "[2/4] Subsetting target loci and line samples..."

# Define temporary intermediate subset VCF path
SUBSET_TEMP_VCF="$TEMP_DIR/intermediate_subset.vcf.gz"

# Check if any filtering options were specified
if [ ${#FILTER_OPTS[@]} -gt 0 ]; then
    
    # Print message executing subset filter
    echo "  -> Applying subset filters with bcftools view..."
    
    # Execute bcftools view with specified filter arguments
    bcftools view "${FILTER_OPTS[@]}" -O z -o "$SUBSET_TEMP_VCF" "$INPUT_VCF"
    
# Otherwise, if no filters were supplied
else
    
    # Print message indicating direct pass-through
    echo "  -> No target or line filters specified. Proceeding with full VCF..."
    
    # Copy input to intermediate subset location
    cp "$INPUT_VCF" "$SUBSET_TEMP_VCF"
    
# End of filter options check
fi

# Sanity Check: Ensure intermediate subset VCF exists and is non-empty
if [ ! -s "$SUBSET_TEMP_VCF" ]; then
    
    # Print error message
    echo "Error: Subsetting step produced an empty intermediate file."
    
    # Exit with error code
    exit 1
    
# End of intermediate file check
fi

# If a target manifest was supplied, verify coordinate + ID pairing (and enforce biallelic REF/ALT when -b is active)
if [ -n "$TARGETS_FILE" ]; then
    
    # Define path for Python-verified/standardized intermediate VCF
    STANDARDIZED_TEMP_VCF="$TEMP_DIR/standardized_subset.vcf"
    
    # Check whether strict canonical biallelic standardization (-b) was requested
    if [ "$BIALLELIC_TARGETS" = true ]; then
        
        # Derive output prefix by stripping .vcf.gz or .vcf extension for the companion QC TSV
        OUTPUT_PREFIX="${OUTPUT_VCF%.vcf.gz}"
        OUTPUT_PREFIX="${OUTPUT_PREFIX%.vcf}"
        QC_REPORT_FILE="${OUTPUT_PREFIX}_biallelic_qc.tsv"
        
        # Print status message for strict biallelic standardization
        echo "  -> Enforcing canonical biallelic REF/ALT alleles and masking off-target calls..."
        
        # Execute standardize_biallelic_targets.py in strict biallelic mode
        "$PYTHON_CMD" "$PY_BIALLELIC_SCRIPT" \
            --input "$SUBSET_TEMP_VCF" \
            --output "$STANDARDIZED_TEMP_VCF" \
            --targets "$TARGETS_FILE" \
            --biallelic \
            --qc-output "$QC_REPORT_FILE"
        
    else
        
        # Print status message for non-destructive coordinate + ID verification
        echo "  -> Verifying exact coordinate + marker ID pairing against target manifest..."
        
        # Execute standardize_biallelic_targets.py in non-destructive coordinate + ID filter mode
        "$PYTHON_CMD" "$PY_BIALLELIC_SCRIPT" \
            --input "$SUBSET_TEMP_VCF" \
            --output "$STANDARDIZED_TEMP_VCF" \
            --targets "$TARGETS_FILE"
        
    # End of biallelic vs non-destructive mode check
    fi
    
    # Compress the verified/standardized VCF back into SUBSET_TEMP_VCF for downstream sorting
    bgzip -c "$STANDARDIZED_TEMP_VCF" > "$SUBSET_TEMP_VCF"
    rm -f "$STANDARDIZED_TEMP_VCF"
    
# End of target manifest Python verification/standardization step
fi

# Calculate count of variant records retained in intermediate subset
SUBSET_VARIANTS=$(bcftools view -H "$SUBSET_TEMP_VCF" 2>/dev/null | wc -l | tr -d ' ' || echo "0")

# If zero variants were retained, issue an informative warning
if [ "$SUBSET_VARIANTS" -eq 0 ]; then
    
    # Print warning message
    echo "  -> Warning: Subsetting retained 0 variant records."
    echo "     None of the target manifest loci were found in '$INPUT_VCF'."
    
else
    
    # Print count of retained records
    echo "  -> Intermediate subset retained $SUBSET_VARIANTS variant records."
    
# End of subset variants check
fi

# Print section 3 header
echo "[3/4] Sorting records in physical coordinate order for TASSEL GUI..."

# Extract non-contig metadata header lines
bcftools view -h "$SUBSET_TEMP_VCF" | grep -v '^##contig=' | grep -v '^#CHROM' > "$TEMP_DIR/new_header.txt"

# Extract canonical chromosomes (Chr1..Chr7) sorted in natural order first
bcftools view -h "$SUBSET_TEMP_VCF" | grep '^##contig=' | grep -E 'ID=(chr|Chr|CHR)[0-9]' | LC_ALL=C sort >> "$TEMP_DIR/new_header.txt" || true

# Extract remaining non-native contigs sorted alphabetically
bcftools view -h "$SUBSET_TEMP_VCF" | grep '^##contig=' | grep -v -E 'ID=(chr|Chr|CHR)[0-9]' | LC_ALL=C sort >> "$TEMP_DIR/new_header.txt" || true

# Append the column header line
bcftools view -h "$SUBSET_TEMP_VCF" | grep '^#CHROM' >> "$TEMP_DIR/new_header.txt"

# Apply the ordered header
bcftools reheader -h "$TEMP_DIR/new_header.txt" -o "$TEMP_DIR/reordered_subset.vcf.gz" "$SUBSET_TEMP_VCF"

# Define temporary sorted VCF path
SORTED_TEMP_VCF="$TEMP_DIR/intermediate_sorted.vcf.gz"

# Execute bcftools sort with specified memory buffer
bcftools sort -m "$SORT_MEMORY" -T "$TEMP_DIR" -O z -o "$SORTED_TEMP_VCF" "$TEMP_DIR/reordered_subset.vcf.gz"

# Sanity Check: Ensure sorted intermediate VCF exists and is non-empty
if [ ! -s "$SORTED_TEMP_VCF" ]; then
    
    # Print error message
    echo "Error: Sorting step produced an empty intermediate file."
    
    # Exit with error code
    exit 1
    
# End of sorted file check
fi

# Print section 4 header
echo "[4/4] Finalizing output VCF and creating genomic indexes..."

# Check if output destination ends in .gz
if [[ "$OUTPUT_VCF" == *.gz ]]; then
    
    # Remove old output files if present
    rm -f "$OUTPUT_VCF" "$OUTPUT_VCF.csi" "$OUTPUT_VCF.tbi" 2>/dev/null || true
    
    # Copy final compressed sorted VCF to destination
    cp -f "$SORTED_TEMP_VCF" "$OUTPUT_VCF"
    
    # Generate CSI index for large chromosome coordinates
    bcftools index -f -c "$OUTPUT_VCF"
    
    # Generate Tabix TBI index for broad viewer compatibility
    bcftools index -f -t "$OUTPUT_VCF" 2>/dev/null || true
    
# Otherwise, if uncompressed output was requested
else
    
    # Remove old uncompressed file if present
    rm -f "$OUTPUT_VCF" 2>/dev/null || true
    
    # Decompress to plain text VCF at destination
    gzip -dc "$SORTED_TEMP_VCF" > "$OUTPUT_VCF"
    
# End of compression format check
fi

# Sanity Check: Ensure final output VCF was created and is non-empty
if [ ! -s "$OUTPUT_VCF" ]; then
    
    # Print error message
    echo "Error: Final output file '$OUTPUT_VCF' was not created or is empty."
    
    # Exit with error code
    exit 1
    
# End of final output file check
fi

# Calculate final site and sample counts for verification
TOTAL_SITES=$(bcftools view -H "$OUTPUT_VCF" 2>/dev/null | wc -l | tr -d ' ' || echo "N/A")
TOTAL_SAMPLES=$(bcftools query -l "$OUTPUT_VCF" 2>/dev/null | wc -l | tr -d ' ' || echo "N/A")

# Print the final completion summary
echo ""
echo "================================================================================"
echo " Pipeline Complete!"
echo " Final Output File:  $OUTPUT_VCF"
if [ -f "$OUTPUT_VCF.csi" ]; then echo " CSI Index:          $OUTPUT_VCF.csi"; fi
if [ -f "$OUTPUT_VCF.tbi" ]; then echo " Tabix Index:        $OUTPUT_VCF.tbi"; fi
if [ "$BIALLELIC_TARGETS" = true ] && [ -f "$QC_REPORT_FILE" ]; then echo " Biallelic QC File:  $QC_REPORT_FILE"; fi
echo " Total Output Sites: $TOTAL_SITES"
echo " Total Output Lines: $TOTAL_SAMPLES"
echo "================================================================================"
