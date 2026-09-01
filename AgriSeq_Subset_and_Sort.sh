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
    echo "  -o, --output         (Required) Path for final output VCF file (.vcf.gz)."
    
    # Print the targets list flag description
    echo "  -t, --targets        (Optional) Target manifest (#CHROM POS ID) or marker list (e.g. V2Targets_manifest.txt)."
    
    # Print the lines list flag description
    echo "  -l, --lines          (Optional) Single-column text file of line/taxa names (e.g. example_lines_list.txt)."
    
    # Print the memory flag description
    echo "  -m, --memory         (Optional) Memory limit for bcftools sorting buffer (default: 16g)."
    
    # Print a blank line
    echo ""
    
    # Print the DETAILS header
    echo "PIPELINE DETAILS:"
    
    # Print description of non-destructive sorting
    echo "  1. Subsets target loci using the exact chromosome coordinates (#CHROM, POS)."
    
    # Print description of sample filtering
    echo "  2. Subsets line/taxa samples using bcftools without dropping sample metadata."
    
    # Print description of TASSEL compatibility
    echo "  3. Sorts records in strict genomic coordinate order so they open cleanly in TASSEL GUI."
    
    # Print description of preservation
    echo "  4. Preserves 100% of headers, contig names (e.g. Chr1A), reference genome REF alleles,"
    
    # Print description of format preservation
    echo "     ALT alleles, and deep sequencing FORMAT/INFO annotations."
    
    # Print a blank line
    echo ""
    
    # Print the REQUIREMENTS header
    echo "REQUIREMENTS:"
    
    # Print the conda requirement description
    echo "  Must be run within a conda environment containing 'bcftools' and 'htslib'."
    
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
    
# End of input file check
fi

# Sanity Check 3: If a targets file is provided, ensure it exists
if [ -n "$TARGETS_FILE" ] && [ ! -f "$TARGETS_FILE" ]; then
    
    # Print an error message
    echo "Error: Targets file '$TARGETS_FILE' does not exist."
    
    # Exit the script with an error code
    exit 1
    
# End of targets file check
fi

# Sanity Check 4: If a lines file is provided, ensure it exists
if [ -n "$LINES_FILE" ] && [ ! -f "$LINES_FILE" ]; then
    
    # Print an error message
    echo "Error: Lines file '$LINES_FILE' does not exist."
    
    # Exit the script with an error code
    exit 1
    
# End of lines file check
fi

# Sanity Check 5: Ensure bcftools is installed and available in PATH
if ! command -v bcftools &> /dev/null; then
    
    # Print an error message
    echo "Error: 'bcftools' is not available. Please activate your conda environment."
    
    # Exit the script with an error code
    exit 1
    
# End of bcftools check
fi

# Sanity Check 6: Ensure bgzip is installed and available in PATH
if ! command -v bgzip &> /dev/null; then
    
    # Print an error message
    echo "Error: 'bgzip' is not available. Please activate your conda environment."
    
    # Exit the script with an error code
    exit 1
    
# End of bgzip check
fi

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
echo "[1/4] Preparing filter parameters and validating files..."

# Initialize array for bcftools filter options
FILTER_OPTS=()

# Check if a targets file was provided
if [ -n "$TARGETS_FILE" ]; then
    
    # Clean targets file into a temporary list stripping Windows carriage returns and empty lines
    grep -ve '^\s*$' "$TARGETS_FILE" | tr -d '\r' > "$TEMP_DIR/cleaned_targets.txt"
    
    # Determine the number of columns in the first non-comment line of the targets file
    FIRST_LINE=$(grep -v '^#' "$TEMP_DIR/cleaned_targets.txt" | head -n 1)
    
    # Count columns in the first data line
    NUM_COLS=$(echo "$FIRST_LINE" | awk '{print NF}')
    
    # Check if targets file has 2 or more columns (e.g. #CHROM POS ID)
    if [ "$NUM_COLS" -ge 2 ]; then
        
        # Count target coordinate loci
        TARGET_COUNT=$(grep -v '^#' "$TEMP_DIR/cleaned_targets.txt" | wc -l | tr -d ' ')
        
        # Print informational message
        echo "  -> Target manifest detected ($NUM_COLS columns, $TARGET_COUNT coordinate loci)."
        
        # Extract the first two columns (CHROM and POS) as tab-delimited for bcftools -T
        awk -v OFS='\t' '!/^#/ {print $1, $2}' "$TEMP_DIR/cleaned_targets.txt" > "$TEMP_DIR/target_regions.tsv"
        
        # Add target coordinate filter flag (-T) and exact start position matching (--targets-overlap 0)
        FILTER_OPTS+=(-T "$TEMP_DIR/target_regions.tsv" --targets-overlap 0)
        
    # Otherwise, if targets file is a 1-column list of marker IDs
    else
        
        # Count target marker IDs
        TARGET_COUNT=$(wc -l < "$TEMP_DIR/cleaned_targets.txt" | tr -d ' ')
        
        # Print informational message
        echo "  -> Single-column marker ID list detected ($TARGET_COUNT marker IDs)."
        
        # Add ID filter flag (-i) to bcftools arguments
        FILTER_OPTS+=(-i 'ID=@'"$TEMP_DIR/cleaned_targets.txt")
        
    # End of target format check
    fi
    
# End of targets file check
fi

# Check if a lines file was provided
if [ -n "$LINES_FILE" ]; then
    
    # Clean line list into a temporary file stripping Windows carriage returns and empty lines
    grep -ve '^\s*$' "$LINES_FILE" | tr -d '\r' > "$TEMP_DIR/cleaned_lines.txt"
    
    # Count the number of lines requested
    REQUESTED_LINES=$(wc -l < "$TEMP_DIR/cleaned_lines.txt" | tr -d ' ')
    
    # Extract all sample names currently present in the input VCF
    bcftools query -l "$INPUT_VCF" > "$TEMP_DIR/all_vcf_samples.txt"
    
    # Intersect requested lines with samples actually present in the input VCF
    awk 'FNR==NR { vcf[$1]=1; next } { if ($1 in vcf) print $1 }' "$TEMP_DIR/all_vcf_samples.txt" "$TEMP_DIR/cleaned_lines.txt" > "$TEMP_DIR/matched_samples.txt"
    
    # Count how many requested lines matched
    MATCHED_COUNT=$(wc -l < "$TEMP_DIR/matched_samples.txt" | tr -d ' ')
    
    # Print status message
    echo "  -> Subsetting $MATCHED_COUNT of $REQUESTED_LINES requested lines present in VCF..."
    
    # Check if any lines matched
    if [ "$MATCHED_COUNT" -eq 0 ]; then
        
        # Print error message
        echo "Error: None of the lines in '$LINES_FILE' matched the sample IDs in '$INPUT_VCF'."
        
        # Print guidance hint
        echo "Please verify that the line designations match the sample names in the VCF header."
        
        # Exit with error code
        exit 1
        
    # End of matched count check
    fi
    
    # Add sample filtering flag to bcftools arguments
    FILTER_OPTS+=(-S "$TEMP_DIR/matched_samples.txt" --force-samples)
    
# End of lines file check
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

# Print section 4 header
echo "[4/4] Finalizing output VCF and creating genomic indexes..."

# Ensure the parent directory for the output file exists
mkdir -p "$(dirname "$OUTPUT_VCF")"

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
