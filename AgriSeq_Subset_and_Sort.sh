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
    echo "  ./AgriSeq_Subset_and_Sort.sh -i <input_vcf> -o <output_vcf> [-t <targets_file>] [-l <lines_file>] [-m <memory>]"
    
    # Print a blank line
    echo ""
    
    # Print the OPTIONS section header
    echo "OPTIONS:"
    
    # Print the help flag description
    echo "  -h, --help           Show this help message and exit."
    
    # Print the input VCF flag description
    echo "  -i, --input          (Required) Path to input VCF file (.vcf or .vcf.gz)."
    
    # Print the output VCF flag description
    echo "  -o, --output         (Required) Path for final output VCF file (.vcf.gz or .vcf)."
    
    # Print the targets list flag description
    echo "  -t, --targets        (Optional) Target manifest with EXACTLY 3 columns (#CHROM POS ID)."
    echo "                       (e.g. V2Targets_manifest.txt or example_target_manifest.txt)."
    echo "                       Single-column ID lists are disallowed due to amplicon indels."
    
    # Print the lines list flag description
    echo "  -l, --lines          (Optional) Single-column text file of line/taxa names (e.g. example_lines_list.txt)."
    
    # Print the memory flag description
    echo "  -m, --memory         (Optional) Memory limit for bcftools sorting buffer (default: 16g)."
    
    # Print a blank line
    echo ""
    
    # Print the DETAILS header
    echo "PIPELINE DETAILS & INTEGRITY CHECKS:"
    
    # Print description of strict manifest validation
    echo "  1. Strict 3-Column Manifest Validation: Enforces exact 3-column (#CHROM POS ID)"
    echo "     format, verifies positive integer coordinates, flags duplicate positions,"
    echo "     and validates target contigs against the reference VCF header. Single-column"
    echo "     marker ID filtering is prohibited because Thermo Fisher TVC VCFs often report"
    echo "     adjacent indels within the amplicon using the same marker ID as the target SNP."
    
    # Print description of sample filtering
    echo "  2. Sample Line Validation: Validates single-column line files, collapses duplicate"
    echo "     taxa, and intersects sample names with the VCF header, reporting missing samples."
    
    # Print description of TASSEL compatibility
    echo "  3. Genomic Ordering: Sorts records in strict physical coordinate order with"
    echo "     natural contig sequencing so VCFs open cleanly in the TASSEL GUI."
    
    # Print description of preservation
    echo "  4. Complete Metadata Preservation: Preserves 100% of headers, contig names"
    echo "     (e.g. Chr1A, TEL10_7E...), REF alleles, ALT alleles, and FORMAT/INFO tags."
    
    # Print a blank line
    echo ""
    
    # Print the REQUIREMENTS header
    echo "REQUIREMENTS:"
    
    # Print the environment requirement description
    echo "  Must be run within an environment containing 'bcftools' and 'htslib' (via Conda or module load)."
    
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

# Sanity Check 6: If a targets file is provided, ensure it exists and is non-empty
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
    
    # Clean targets file into a temporary list stripping Windows carriage returns and empty lines
    grep -ve '^\s*$' "$TARGETS_FILE" | tr -d '\r' > "$TEMP_DIR/cleaned_targets.txt"
    
    # Count total non-comment data lines in targets file
    TARGET_COUNT=$(grep -v '^#' "$TEMP_DIR/cleaned_targets.txt" | wc -l | tr -d ' ')
    
    # Sanity Check: Ensure target manifest has at least one valid data row
    if [ "$TARGET_COUNT" -eq 0 ]; then
        
        # Print error message
        echo "Error: Target manifest '$TARGETS_FILE' contains no data records (only comments or whitespace)."
        
        # Print guidance hint
        echo "Please provide a valid target manifest containing '#CHROM POS ID' rows."
        
        # Exit with error code
        exit 1
        
    # End of target count check
    fi
    
    # Sanity Check: STRICT 3-COLUMN CHECK
    # Check if any non-comment line has a column count different from 3
    INVALID_COL_LINE=$(awk '!/^#/ {if (NF != 3) {print NR, NF, $0; exit}}' "$TEMP_DIR/cleaned_targets.txt")
    
    # If a line with a column count other than 3 was detected
    if [ -n "$INVALID_COL_LINE" ]; then
        
        # Extract the line number of the offending row
        BAD_LINE_NUM=$(echo "$INVALID_COL_LINE" | awk '{print $1}')
        
        # Extract the number of columns found on that row
        FOUND_COLS=$(echo "$INVALID_COL_LINE" | awk '{print $2}')
        
        # Extract the content of the offending row
        BAD_LINE_TEXT=$(echo "$INVALID_COL_LINE" | cut -d' ' -f3-)
        
        # Print descriptive error banner
        echo "================================================================================"
        echo "Error: Target manifest format violation in '$TARGETS_FILE'!"
        echo "Found $FOUND_COLS columns at line $BAD_LINE_NUM: '$BAD_LINE_TEXT'"
        echo ""
        echo "The target manifest MUST contain EXACTLY 3 whitespace-separated columns:"
        echo "  Column 1: #CHROM (Contig or chromosome name, e.g. Chr1A or TEL10_7E739933655)"
        echo "  Column 2: POS    (1-based integer physical position, e.g. 12345678)"
        echo "  Column 3: ID     (Target marker identifier, e.g. AX-12345678)"
        echo ""
        echo "WHY EXACT 3 COLUMNS ARE REQUIRED:"
        echo "Passing a single-column marker name list is strictly prohibited because"
        echo "Thermo Fisher TVC VCFs frequently report multiple off-target amplicon indels"
        echo "at the same locus that share the exact same marker ID as the target SNP."
        echo "Filtering strictly by coordinates (#CHROM, POS) is mandatory to isolate"
        echo "the true target variant and reject adjacent amplicon indels."
        echo ""
        echo "Example valid manifest lines:"
        echo "  #CHROM             POS   ID"
        echo "  Chr1A              1169383  CHS21_1A001169383"
        echo "  TEL10_7E739933655  255      TEL10_7E739933655"
        echo "================================================================================"
        
        # Exit with error code
        exit 1
        
    # End of strict 3-column check
    fi
    
    # Sanity Check: POS column must be a positive integer on every data row
    NON_NUMERIC_POS=$(awk '!/^#/ {if ($2 !~ /^[0-9]+$/ || $2 <= 0) {print NR, $2, $0; exit}}' "$TEMP_DIR/cleaned_targets.txt")
    
    # If a non-numeric or non-positive position was detected
    if [ -n "$NON_NUMERIC_POS" ]; then
        
        # Extract the line number
        BAD_POS_LINE=$(echo "$NON_NUMERIC_POS" | awk '{print $1}')
        
        # Extract the invalid position value
        BAD_POS_VAL=$(echo "$NON_NUMERIC_POS" | awk '{print $2}')
        
        # Print error message
        echo "Error: Invalid physical position '$BAD_POS_VAL' at line $BAD_POS_LINE in '$TARGETS_FILE'."
        echo "Column 2 (POS) must be a positive integer genomic coordinate."
        
        # Exit with error code
        exit 1
        
    # End of POS numeric check
    fi
    
    # Sanity Check: Detect duplicate (CHROM, POS) coordinates in the manifest
    DUP_COORDS=$(awk '!/^#/ {print $1, $2}' "$TEMP_DIR/cleaned_targets.txt" | LC_ALL=C sort | uniq -c | awk '$1 > 1 {print $2 ":" $3 " (" $1 " entries)"}')
    
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
        
        # Explain handling of duplicate coordinates
        echo "     Bcftools will subset each unique coordinate once."
        
    # End of duplicate coordinates check
    fi
    
    # Sanity Check: Verify target manifest chromosomes exist in the input VCF header
    awk '!/^#/ {print $1}' "$TEMP_DIR/cleaned_targets.txt" | LC_ALL=C sort -u > "$TEMP_DIR/manifest_chroms.txt"
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
    
    # Print status message confirming validated target manifest
    echo "  -> Validated 3-column target manifest: $TARGET_COUNT loci across $TOTAL_MANIFEST_CHROMS contigs."
    
    # Extract the first two columns (CHROM and POS) as tab-delimited for bcftools -T
    awk -v OFS='\t' '!/^#/ {print $1, $2}' "$TEMP_DIR/cleaned_targets.txt" > "$TEMP_DIR/target_regions.tsv"
    
    # Add target coordinate filter flag (-T) and exact start position matching (--targets-overlap 0)
    FILTER_OPTS+=(-T "$TEMP_DIR/target_regions.tsv" --targets-overlap 0)
    
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

# Calculate count of variant records retained in intermediate subset
SUBSET_VARIANTS=$(bcftools view -H "$SUBSET_TEMP_VCF" 2>/dev/null | wc -l | tr -d ' ' || echo "0")

# If zero variants were retained, issue an informative warning
if [ "$SUBSET_VARIANTS" -eq 0 ]; then
    
    # Print warning message
    echo "  -> Warning: Subsetting retained 0 variant records."
    echo "     None of the target manifest coordinates were found in '$INPUT_VCF'."
    
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
echo " Total Output Sites: $TOTAL_SITES"
echo " Total Output Lines: $TOTAL_SAMPLES"
echo "================================================================================"
