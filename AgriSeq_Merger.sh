#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define a function to display the help and usage section
function show_help {
    
    # Print the header line
    echo "================================================================================"
    
    # Print the script title
    echo " AgriSeq_Merger.sh - A robust pipeline for processing and merging AgriSeq VCFs"
    
    # Print the footer line
    echo "================================================================================"
    
    # Print a blank line
    echo ""
    
    # Print the USAGE section header
    echo "USAGE:"
    
    # Print the general command structure
    echo "  ./AgriSeq_Merger.sh -z <zip_dir> -r <ref_index> -c <contigs> -o <output_vcf> [-b <batch_size>] [-m <map_file>] [-p <non_native_prefix>]"
    
    # Print a blank line
    echo ""
    
    # Print the OPTIONS section header
    echo "OPTIONS:"
    
    # Print the help flag description
    echo "  -h, --help           Show this help message and exit."
    
    # Print the zip directory flag description
    echo "  -z, --zip-dir        Directory containing the original .zip files."
    
    # Print the reference index flag description
    echo "  -r, --ref-index      Reference index string to prefix markers (e.g., CHS21)."
    
    # Print the contigs flag description
    echo "  -c, --contigs        Comma-separated list of contigs belonging to the reference."
    
    # Print a continuation of the contigs description
    echo "                       (e.g., Chr1A,Chr1B,Chr1D,...)"
    
    # Print the output file flag description
    echo "  -o, --output         Name of the final merged output file (e.g., Final_Merged.vcf.gz)."
    
    # Print the batch size flag description
    echo "  -b, --batch-size     (Optional) Number of VCFs to merge simultaneously (default: 100)."
    
    # Print the map file flag description
    echo "  -m, --map-file       (Optional) Text file mapping sample IDs to line names (old_name new_name) to reheader the final VCF."
    
    # Print the non-native prefix flag description
    echo "  -p, --non-native     (Optional) String to concatenate prior to position on non-native contigs."
    
    # Print a blank line
    echo ""
    
    # Print the Python script usage header
    echo "PYTHON SCRIPT USAGE (rename_vcf_markers.py):"
    
    # Print the explanation that the python script is called internally
    echo "  The pipeline internally calls rename_vcf_markers.py. You can also run it"
    
    # Print the explanation for viewing python help
    echo "  independently. View its help section by running: python rename_vcf_markers.py -h"
    
    # Print a blank line
    echo ""
    
    # Print the REQUIREMENTS header
    echo "REQUIREMENTS:"
    
    # Print the conda requirement description
    echo "  Must be run within a conda environment containing 'bcftools', 'htslib', and 'python'."
    
    # Print the final border line
    echo "================================================================================"
    
# End of the show_help function
}

# Initialize the ZIP_DIR variable
ZIP_DIR=""

# Initialize the REF_INDEX variable
REF_INDEX=""

# Initialize the CONTIGS variable
CONTIGS=""

# Initialize the MERGED_OUTPUT variable
MERGED_OUTPUT=""

# Initialize the BATCH_SIZE variable with a default of 100
BATCH_SIZE=100

# Initialize the MAP_FILE variable
MAP_FILE=""

# Initialize the NON_NATIVE_PREFIX variable
NON_NATIVE_PREFIX=""

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
        # If the zip directory flag is passed
        -z|--zip-dir)
            
            # Assign the next argument to ZIP_DIR
            ZIP_DIR="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of zip directory case
            ;;
            
        # If the reference index flag is passed
        -r|--ref-index)
            
            # Assign the next argument to REF_INDEX
            REF_INDEX="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of reference index case
            ;;
            
        # If the contigs list flag is passed
        -c|--contigs)
            
            # Assign the next argument to CONTIGS
            CONTIGS="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of contigs list case
            ;;
            
        # If the output file flag is passed
        -o|--output)
            
            # Assign the next argument to MERGED_OUTPUT
            MERGED_OUTPUT="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of output file case
            ;;
            
        # If the batch size flag is passed
        -b|--batch-size)
            
            # Assign the next argument to BATCH_SIZE
            BATCH_SIZE="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of batch size case
            ;;
            
        # If the map file flag is passed
        -m|--map-file)
            
            # Assign the next argument to MAP_FILE
            MAP_FILE="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of map file case
            ;;
            
        # If the non-native prefix flag is passed
        -p|--non-native)
            
            # Assign the next argument to NON_NATIVE_PREFIX
            NON_NATIVE_PREFIX="$2"
            
            # Shift the arguments by 2
            shift 2
            
            # End of non native prefix case
            ;;
            
        # If the deduplicate ids flag is passed
        -d|--deduplicate-ids)
            
            # Enable deduplication logic
            DEDUPLICATE_IDS=true
            
            # Shift the arguments by 1
            shift 1
            
            # End of deduplicate ids case
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

# Sanity Check 1: Check if any of the required arguments are empty
if [ -z "$ZIP_DIR" ] || [ -z "$REF_INDEX" ] || [ -z "$CONTIGS" ] || [ -z "$MERGED_OUTPUT" ]; then
    
    # Print an error indicating missing arguments
    echo "Error: Missing required arguments."
    
    # Show the help section to guide the user
    show_help
    
    # Exit the script with an error code
    exit 1
    
# End of the required argument check
fi

# Sanity Check 2: Ensure the provided zip directory actually exists
if [ ! -d "$ZIP_DIR" ]; then
    
    # Print an error message
    echo "Error: Zip directory '$ZIP_DIR' does not exist or is not a directory."
    
    # Exit the script with an error code
    exit 1
    
# End of zip directory check
fi

# Sanity Check 3: Ensure python is installed and available in PATH
if ! command -v python &> /dev/null; then
    
    # Print an error message
    echo "Error: 'python' is not installed or not available in your PATH."
    
    # Exit the script with an error code
    exit 1
    
# End of python check
fi

# Sanity Check 4: Ensure bcftools is installed and available in PATH
if ! command -v bcftools &> /dev/null; then
    
    # Print an error message
    echo "Error: 'bcftools' is not available. Please activate your conda environment."
    
    # Exit the script with an error code
    exit 1
    
# End of bcftools check
fi

# Sanity Check 5: Ensure bgzip is installed and available in PATH
if ! command -v bgzip &> /dev/null; then
    
    # Print an error message
    echo "Error: 'bgzip' is not available. Please activate your conda environment."
    
    # Exit the script with an error code
    exit 1
    
# End of bgzip check
fi

# Sanity Check 6: If a map file is provided, ensure it actually exists
if [ -n "$MAP_FILE" ] && [ ! -f "$MAP_FILE" ]; then
    
    # Print an error message
    echo "Error: Map file '$MAP_FILE' does not exist or is not a valid file."
    
    # Exit the script with an error code
    exit 1
    
# End of map file check
fi

# Create a secure temporary directory managed by the OS
TEMP_DIR=$(mktemp -d)

# Set a trap to ensure the temporary directory is safely deleted whenever the script exits
trap 'rm -rf "$TEMP_DIR"' EXIT

# Print a message for the setup step
echo "[1/5] Setting up directories..."

# Define the output directory name for processed VCFs inside the temp directory
FINAL_VCF_DIR="$TEMP_DIR/Processed_Ready_VCFs"

# Create the final processed directory
mkdir -p "$FINAL_VCF_DIR"

# Define the temporary directory name for unzipping files inside our trapped temp directory
TEMP_UNZIP_DIR="$TEMP_DIR/temp_unzipped_vcfs"

# Create the temporary unzipped directory
mkdir -p "$TEMP_UNZIP_DIR"

# Print a message for the extraction step
echo "[2/5] Extracting original .zip files from $ZIP_DIR..."

# Find all zip files and loop through them
find "$ZIP_DIR" -type f -name "*.zip" | while read -r zip_file; do
    
    # Extract the base name of the zip file
    zip_base=$(basename "$zip_file" .zip)
    
    # Create a unique subdirectory for this zip file inside the temp dir
    mkdir -p "$TEMP_UNZIP_DIR/$zip_base"
    
    # Extract the zip file quietly into its unique subdirectory
    unzip -q -o "$zip_file" -d "$TEMP_UNZIP_DIR/$zip_base"
    
# End of the find and unzip loop
done

# Print a message for the processing step
echo "[3/5] Renaming markers, bgzipping, and indexing..."

# Find all extracted vcf.gz files and loop through them
find "$TEMP_UNZIP_DIR" -type f -name "*.vcf.gz" | while read -r vcf_file; do
    
    # Extract the base filename by stripping the path and the .gz extension
    vcf_base=$(basename "$vcf_file" .gz)
    
    # Extract the parent directory name (which is the unique zip_base we created)
    parent_dir=$(basename "$(dirname "$vcf_file")")
    
    # Create a globally unique base name by combining the zip name and the vcf name
    unique_base="${parent_dir}_${vcf_base}"
    
    # Initialize the python command as a string using the globally unique output name
    PY_CMD="python rename_vcf_markers.py --input \"$vcf_file\" --output \"$FINAL_VCF_DIR/$unique_base\" --ref-index \"$REF_INDEX\" --ref-contigs \"$CONTIGS\""
    
    # Check if a non-native prefix was provided
    if [ -n "$NON_NATIVE_PREFIX" ]; then
        
        # Append the non-native prefix flag to the python command
        PY_CMD="$PY_CMD --non-native-prefix \"$NON_NATIVE_PREFIX\""
        
    # End of the non-native prefix check
    fi
    
    # Run the constructed Python command
    eval $PY_CMD
    
    # Compress the newly created VCF using bgzip
    bgzip -f "$FINAL_VCF_DIR/$unique_base"
    
    # Create a CSI index for the bgzipped file to support large chromosomes
    bcftools index -f -c "$FINAL_VCF_DIR/$unique_base.gz"
    
# End of the VCF processing loop
done

# Print a message for the batch merging step
echo "[4/5] Batch merging VCFs to avoid open-file limits..."

# Define a temporary directory name for intermediate batch files inside our trapped temp directory
BATCH_DIR="$TEMP_DIR/temp_batch_merges"

# Create the batch directory
mkdir -p "$BATCH_DIR"

# Load all processed VCF file paths into an array
VCF_FILES=("$FINAL_VCF_DIR"/*.vcf.gz)

# Count the total number of files in the array
TOTAL_FILES=${#VCF_FILES[@]}

# Calculate the total number of batches
TOTAL_BATCHES=$(( (TOTAL_FILES + BATCH_SIZE - 1) / BATCH_SIZE ))

# Print the number of files discovered, the batch size, and total batches
echo "Found $TOTAL_FILES files. Merging in $TOTAL_BATCHES batches of $BATCH_SIZE..."

# Initialize a counter for the current batch
BATCH_COUNT=0

# Loop through the files in increments of the batch size
for ((i=0; i < TOTAL_FILES; i+=BATCH_SIZE)); do
    
    # Increment the batch counter
    BATCH_COUNT=$((BATCH_COUNT + 1))
    
    # Define the path to the text file that will hold the list of VCFs for this batch
    BATCH_LIST="$BATCH_DIR/batch_${BATCH_COUNT}_list.txt"
    
    # Define the output VCF name for this specific batch
    BATCH_OUTPUT="$BATCH_DIR/batch_${BATCH_COUNT}.vcf.gz"
    
    # Loop through the specific files that belong to this batch
    for ((j=i; j < i+BATCH_SIZE && j < TOTAL_FILES; j++)); do
        
        # Append the current file path to the batch list
        echo "${VCF_FILES[j]}" >> "$BATCH_LIST"
        
    # End of the inner loop for generating the batch list
    done
    
    # Print the current batch being merged out of the total batches
    echo "  -> Merging batch $BATCH_COUNT of $TOTAL_BATCHES..."
    
    # Run bcftools merge on the list of files to generate the batch output
    bcftools merge -l "$BATCH_LIST" -O z -o "$BATCH_OUTPUT"
    
    # Create a CSI index for the newly created intermediate batch file
    bcftools index -f -c "$BATCH_OUTPUT"
    
# End of the outer loop for processing all batches
done

# Check if exactly one batch was needed
if [ "$BATCH_COUNT" -eq 1 ]; then
    
    # Print a message indicating a direct rename
    echo "Only 1 batch was needed. Renaming to final output..."
    
    # Move the single batch file to the final merged output location
    mv "$BATCH_DIR/batch_1.vcf.gz" "$MERGED_OUTPUT"
    
    # Move the single batch index file to the final merged output location
    mv "$BATCH_DIR/batch_1.vcf.gz.csi" "$MERGED_OUTPUT.csi"
    
# Otherwise, if multiple batches were created
else
    
    # Print a message indicating the final master merge
    echo "Merging all $BATCH_COUNT intermediate batches into $MERGED_OUTPUT..."
    
    # Find all intermediate batch VCFs and write them to a final text list
    find "$BATCH_DIR" -type f -name "*.vcf.gz" > "$TEMP_DIR/final_batch_list.txt"
    
    # Run bcftools merge on the intermediate batches to create the final unified VCF
    bcftools merge -l "$TEMP_DIR/final_batch_list.txt" -O z -o "$MERGED_OUTPUT"
    
    # Create a CSI index for the final merged VCF
    bcftools index -f -c "$MERGED_OUTPUT"
    
# End of the single vs multiple batch check
fi

# Check if the user requested deduplication of marker IDs
if [ "$DEDUPLICATE_IDS" = true ]; then
    
    # Print a message indicating deduplication is happening
    echo "[*] Deduplicating concatenated marker IDs (keeping first ID only)..."
    
    # Use bcftools view and awk to strip everything after the first semicolon in the ID column
    bcftools view "$MERGED_OUTPUT" | awk -F'\t' -v OFS='\t' '/^#/ {print; next} { split($3, a, ";"); $3 = a[1]; print }' | bgzip -c > "$TEMP_DIR/dedup.vcf.gz"
    
    # Overwrite the merged VCF with the cleaned VCF
    mv "$TEMP_DIR/dedup.vcf.gz" "$MERGED_OUTPUT"
    
    # Re-index the cleaned VCF
    bcftools index -f -c "$MERGED_OUTPUT"
    
# End of the deduplication check
fi

# Check if a mapping file was provided to reheader the VCF
if [ -n "$MAP_FILE" ]; then
    
    # Print a message for the reheadering step
    echo "[5/5] Reheadering sample line names using the map file..."
    
    # Extract the actual valid sample IDs that exist inside the merged VCF
    bcftools query -l "$MERGED_OUTPUT" > "$TEMP_DIR/vcf_samples.txt"
    
    # Use awk to generate a filtered mapping file containing only those samples present in the VCF
    awk 'FNR==NR { vcf[$1]=1; next } { if ($1 in vcf) print $1, $2 }' "$TEMP_DIR/vcf_samples.txt" "$MAP_FILE" > "$TEMP_DIR/valid_mapping.txt"
    
    # Apply the filtered map file to reheader the VCF, outputting to a temporary file
    bcftools reheader -s "$TEMP_DIR/valid_mapping.txt" -o "$TEMP_DIR/reheaded.vcf.gz" "$MERGED_OUTPUT"
    
    # Replace the original merged VCF with the newly reheaded VCF
    mv "$TEMP_DIR/reheaded.vcf.gz" "$MERGED_OUTPUT"
    
    # Re-index the final file with CSI since the content changed
    bcftools index -f -c "$MERGED_OUTPUT"
    
# Otherwise, skip the reheadering step
else
    
    # Print a message indicating reheadering is skipped
    echo "[5/5] No map file provided. Skipping reheader..."
    
# End of the map file check
fi

# Print a final message about the trap cleaning up
echo "Temporary files are being automatically cleaned up by the trap..."

# Print a final success header
echo "=========================================="

# Print a final success statement
echo "Pipeline complete!"

# Print the location of the final merged file
echo "Final merged file: $MERGED_OUTPUT"

# Print a final success footer
echo "=========================================="
