#!/usr/bin/env python3
"""
quick_multiplex.py - Helper module for AgriSeq_Merge_Subset_Sort_Quick.sh

Provides high-speed helper routines for:
  1. index-samples   : Discovering raw sample VCFs from zip archives or directory,
                       mapping barcodes to cultivar names, and filtering lines.
  2. prepare-manifest: Validating and sorting canonical 5-column target manifests
                       into user-defined or standard contig order.
  3. project-sample  : Evaluating a single sample's raw VCF against canonical
                       manifest targets and emitting a uniform single-column file.
  4. assemble-matrix : Horizontally pasting sample columns with manifest metadata.
  5. sort-vcf        : Canonical genomic sorting for TASSEL GUI compatibility.
  6. generate-qc     : Computing per-marker call statistics and companion QC TSV.
"""

# Import standard argparse module for command line arguments
import argparse

# Import glob module for directory path wildcard matching
import glob

# Import gzip module for reading and writing compressed BGZF streams
import gzip

# Import os module for filesystem paths and existence checks
import os

# Import re module for natural alphanumeric string splitting and regex matching
import re

# Import sys module for system streams and exit codes
import sys

# Import zipfile module for reading provider zip archives directly
import zipfile


# Define canonical ordered wheat chromosome list (Chr1A through Chr7D)
WHEAT_ORDER = [f"Chr{g}{s}" for g in range(1, 8) for s in ("A", "B", "D")]

# Build dictionary mapping canonical wheat chromosome name to sort rank index
WHEAT_RANK = {c: i for i, c in enumerate(WHEAT_ORDER)}


def natural_sort_key(s):
    """Sort key for natural alphanumeric string sorting (e.g. Chr2 before Chr10)."""
    
    # Split string into digit and non-digit tokens and cast numeric tokens to integers
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r"(\d+)", s)]


def load_contig_order(contigs_path):
    """Load ordered contig list from a headless text file (one contig per line)."""
    
    # Check if a contigs file path was provided and exists on disk
    if not contigs_path or not os.path.isfile(contigs_path):
        
        # Return empty list if no valid file was supplied
        return []

    # Initialize list to hold ordered contig names
    contigs = []

    # Open headless contigs file for reading
    with open(contigs_path, "r", encoding="utf-8-sig") as f:
        
        # Iterate through lines in the file
        for line in f:
            
            # Strip whitespace and line breaks
            c = line.strip()
            
            # Add non-empty lines that do not start with comments
            if c and not c.startswith("#"):
                
                # Append contig name to list
                contigs.append(c)

    # Return parsed contigs list
    return contigs


def build_contig_rank(custom_contigs=None, seen_contigs=None):
    """Construct contig rank dictionary mapping contig names to sorting priorities."""
    
    # Check if user provided custom contigs
    if custom_contigs:
        
        # Build dictionary from custom user-supplied contig order
        rank = {c: i for i, c in enumerate(custom_contigs)}
        
        # Return custom rank dictionary and list
        return rank, list(custom_contigs)

    # Check if any seen contigs match canonical wheat chromosomes
    if seen_contigs and any(c in WHEAT_RANK for c in seen_contigs):
        
        # Return default wheat rank dictionary and canonical wheat chromosome list
        return dict(WHEAT_RANK), list(WHEAT_ORDER)

    # Default to empty rank for non-wheat species without custom file
    return {}, []


def trim_minimal(p, r, a):
    """Normalize alleles using minimal parsimony trimming."""
    
    # Convert reference and alternate alleles to uppercase
    r, a = r.upper(), a.upper()

    # Trim matching trailing bases while retaining at least one base
    while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
        
        # Remove trailing base from both alleles
        r, a = r[:-1], a[:-1]

    # Trim matching leading bases and increment genomic coordinate
    while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        
        # Remove leading base from both alleles
        r, a = r[1:], a[1:]
        
        # Increment genomic position by one
        p += 1

    # Return normalized position, reference allele, and alternate allele
    return p, r, a


def sort_key_record(chrom, pos, contig_rank=None):
    """Sort key for genomic records using contig rank and natural sorting."""
    
    # Check if contig rank dictionary was provided and contig exists in rank
    if contig_rank and chrom in contig_rank:
        
        # Return tuple with priority 0, contig rank, position, and chromosome
        return (0, contig_rank[chrom], pos, chrom)

    # Return tuple with priority 1, natural sort key of contig, position, and chromosome
    return (1, natural_sort_key(chrom), pos, chrom)


# ==============================================================================
# Subcommand 1: index-samples
# ==============================================================================
def cmd_index_samples(args):
    
    # Initialize dictionary mapping raw sample barcode IDs to line names
    id_map = {}

    # Check if a sample mapping file was provided and exists on disk
    if args.map_file and os.path.exists(args.map_file):
        
        # Open the sample mapping file for reading
        with open(args.map_file, "r") as mf:
            
            # Iterate through each line of the mapping file
            for line in mf:
                
                # Strip leading and trailing whitespace
                s = line.strip()
                
                # Check if the line contains non-empty text
                if s:
                    
                    # Split line into whitespace-delimited columns
                    parts = s.split()
                    
                    # Ensure at least two columns exist (raw_id and line_name)
                    if len(parts) >= 2:
                        
                        # Store mapping from raw ID to designated line name
                        id_map[parts[0]] = parts[1]

    # Initialize set of requested line names if filtering is requested
    target_lines = None

    # Check if a lines subset filter file was provided and exists on disk
    if args.lines_file and os.path.exists(args.lines_file):
        
        # Create an empty set for target line names
        target_lines = set()
        
        # Open the lines list file for reading
        with open(args.lines_file, "r") as lf:
            
            # Iterate through each line of the lines list file
            for line in lf:
                
                # Strip leading and trailing whitespace
                s = line.strip()
                
                # Check if line contains non-empty text
                if s:
                    
                    # Add requested line name to target set
                    target_lines.add(s)

    # Initialize list to hold discovered sample descriptor tuples
    discovered = []

    # Define inner helper function to extract sample barcode ID from VCF header
    def extract_sample_id_from_stream(stream):
        
        # Iterate through header lines of the VCF stream
        for line in stream:
            
            # Check if line is the #CHROM header line
            if line.startswith("#CHROM"):
                
                # Split #CHROM header line into tab-delimited columns
                parts = line.strip().split("\t")
                
                # Ensure the 10th column (sample identifier) is present
                if len(parts) >= 10:
                    
                    # Return sample barcode ID from column 10
                    return parts[9]
                
                # Break loop once #CHROM line has been evaluated
                break
        
        # Return None if no sample ID could be identified
        return None

    # Search for raw single-sample VCFs inside zip archives if zip directory provided
    if args.zip_dir and os.path.isdir(args.zip_dir):
        
        # Iterate through all zip archives in sorted order
        for zf_path in sorted(glob.glob(os.path.join(args.zip_dir, "*.zip"))):
            
            # Wrap archive inspection in exception handling block
            try:
                
                # Open zip archive in read mode
                with zipfile.ZipFile(zf_path, "r") as z:
                    
                    # Iterate through all entries inside current zip archive
                    for entry in sorted(z.namelist()):
                        
                        # Check if entry is a compressed VCF and not an index
                        if entry.endswith(".vcf.gz") and not entry.endswith(".tbi"):
                            
                            # Open inner compressed VCF entry directly from zip
                            with z.open(entry) as raw_file, gzip.open(raw_file, "rt") as gz_file:
                                
                                # Extract sample identifier from VCF header
                                sample_id = extract_sample_id_from_stream(gz_file)
                            
                            # Check if valid sample identifier was extracted
                            if sample_id:
                                
                                # Resolve line name using map file or default to barcode ID
                                line_name = id_map.get(sample_id, sample_id)
                                
                                # Check if sample passes line subset filter criteria
                                if (
                                    target_lines is None
                                    or line_name in target_lines
                                    or sample_id in target_lines
                                ):
                                    
                                    # Append discovered sample descriptor tuple
                                    discovered.append(("ZIP", zf_path, entry, sample_id, line_name))
            
            # Catch exceptions during zip processing and emit warning
            except Exception as e:
                
                # Print warning message to stderr
                sys.stderr.write(f"Warning: Failed reading {zf_path}: {e}\n")

    # Search for pre-extracted raw VCFs if input directory provided
    if args.input_dir and os.path.isdir(args.input_dir):
        
        # Search recursively for all .vcf.gz files in input directory
        for vcf_path in sorted(glob.glob(os.path.join(args.input_dir, "**", "*.vcf.gz"), recursive=True)):
            
            # Skip any index files ending in .tbi
            if vcf_path.endswith(".tbi"):
                
                # Continue to next file
                continue
            
            # Wrap file inspection in exception handling block
            try:
                
                # Open compressed VCF file for text reading
                with gzip.open(vcf_path, "rt") as gz_file:
                    
                    # Extract sample identifier from VCF header
                    sample_id = extract_sample_id_from_stream(gz_file)
                
                # Check if valid sample identifier was extracted
                if sample_id:
                    
                    # Resolve line name using map file or default to barcode ID
                    line_name = id_map.get(sample_id, sample_id)
                    
                    # Check if sample passes line subset filter criteria
                    if (
                        target_lines is None
                        or line_name in target_lines
                        or sample_id in target_lines
                    ):
                        
                        # Append discovered sample descriptor tuple
                        discovered.append(("FILE", vcf_path, "", sample_id, line_name))
            
            # Catch exceptions during file processing and emit warning
            except Exception as e:
                
                # Print warning message to stderr
                sys.stderr.write(f"Warning: Failed reading {vcf_path}: {e}\n")

    # Initialize set to track line names already registered
    seen_lines = set()

    # Initialize list for deduplicated sample descriptor tuples
    unique_discovered = []

    # Iterate through all discovered sample descriptors
    for item in discovered:
        
        # Extract designated line name
        lname = item[4]
        
        # Check if line name has already been registered
        if lname in seen_lines:
            
            # Emit notice to stderr about duplicate entry
            sys.stderr.write(f"Notice: Duplicate entry for line '{lname}' encountered. Keeping first.\n")
            
            # Skip duplicate entry
            continue
        
        # Register line name as seen
        seen_lines.add(lname)
        
        # Add sample descriptor to unique list
        unique_discovered.append(item)

    # Open output index manifest file for writing
    with open(args.output, "w") as out_f:
        
        # Iterate through unique discovered sample descriptors
        for item in unique_discovered:
            
            # Write tab-delimited descriptor line
            out_f.write("\t".join(item) + "\n")

    # Print summary message indicating number of discovered samples
    print(f"  -> Discovered {len(unique_discovered)} matching samples for multiplexing.")

    # Check if line subset filtering was active
    if target_lines is not None:
        
        # Compute set difference of missing requested lines
        missing = target_lines - seen_lines
        
        # Check if any requested lines were absent from input VCFs
        if missing:
            
            # Print notice regarding missing requested lines
            print(f"  -> Notice: {len(missing)} requested lines absent from input VCFs (e.g. {sorted(list(missing))[:5]}...)")


# ==============================================================================
# Subcommand 2: prepare-manifest
# ==============================================================================
def cmd_prepare_manifest(args):
    
    # Load custom contig ordering from file if provided by user
    custom_contigs = load_contig_order(args.contigs_file) if getattr(args, "contigs_file", None) else []

    # Open target manifest file for reading
    with open(args.manifest, "r") as f:
        
        # Read header line and strip newline characters
        hdr_line = f.readline().rstrip("\r\n")
        
        # Strip comment characters and split header by tabs
        cols = [c.strip("# ") for c in hdr_line.split("\t")]
        
        # Convert all header column names to uppercase
        u_cols = [c.upper() for c in cols]

        # Locate column indexes for required fields with exception handling
        try:
            
            # Identify chromosome column index
            c_idx = u_cols.index("CHROM") if "CHROM" in u_cols else u_cols.index("CHR")
            
            # Identify position column index
            p_idx = u_cols.index("POS")
            
            # Identify marker ID column index
            i_idx = u_cols.index("ID")
            
            # Identify reference allele column index
            r_idx = u_cols.index("REF")
            
            # Identify alternate allele column index
            a_idx = u_cols.index("ALT")
            
        # Catch missing column exception
        except ValueError:
            
            # Print error message with column names found
            sys.stderr.write(f"Error: Manifest must contain #CHROM, POS, ID, REF, ALT. Found: {cols}\n")
            
            # Exit script with error code 1
            exit(1)

        # Initialize list to hold target record tuples
        records = []

        # Initialize set to record unique contig names encountered
        seen_contigs = set()

        # Iterate through data rows in the target manifest
        for line in f:
            
            # Strip trailing newline characters
            s = line.rstrip("\r\n")
            
            # Skip empty lines and comment lines
            if not s or s.startswith("#"):
                
                # Continue to next line
                continue
            
            # Split row into tab-delimited columns
            parts = s.split("\t")
            
            # Ensure line contains enough columns for all required indexes
            if len(parts) <= max(c_idx, p_idx, i_idx, r_idx, a_idx):
                
                # Skip truncated line
                continue
            
            # Extract contig name
            chrom = parts[c_idx].strip()
            
            # Extract genomic coordinate as integer
            pos = int(parts[p_idx].strip())
            
            # Extract marker ID
            mid = parts[i_idx].strip()
            
            # Extract reference allele in uppercase
            ref = parts[r_idx].strip().upper()
            
            # Extract alternate allele in uppercase
            alt = parts[a_idx].strip().upper()
            
            # Append target record tuple
            records.append((chrom, pos, mid, ref, alt))
            
            # Add contig name to seen contigs set
            seen_contigs.add(chrom)

    # Build contig ranking dictionary and primary ordered list based on custom file or detection
    contig_rank, primary_order = build_contig_rank(custom_contigs, seen_contigs)

    # Sort target records using resolved contig rank and natural sorting
    sorted_records = sorted(records, key=lambda r: sort_key_record(r[0], r[1], contig_rank))

    # Open output targets file for writing clean canonical loci
    with open(args.output_targets, "w") as out_t:
        
        # Iterate through sorted target records
        for r in sorted_records:
            
            # Write tab-delimited target record line
            out_t.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\n")

    # Order contigs with primary order contigs first (if seen in manifest)
    ordered_contigs = [c for c in primary_order if c in seen_contigs]

    # Collect and naturally sort any extra contigs present in manifest but absent from primary rank
    extra_contigs = sorted([c for c in seen_contigs if c not in contig_rank], key=natural_sort_key)

    # Append extra naturally sorted contigs to ordered contigs list
    ordered_contigs.extend(extra_contigs)

    # Open output contigs file for writing ordered contig definitions
    with open(args.output_contigs, "w") as out_c:
        
        # Iterate through ordered contig names
        for c in ordered_contigs:
            
            # Write contig name line
            out_c.write(c + "\n")

    # Print status message confirming validated target and contig counts
    print(f"  -> Validated {len(sorted_records)} canonical manifest targets across {len(ordered_contigs)} contigs.")


# ==============================================================================
# Subcommand 3: project-sample
# ==============================================================================
def cmd_project_sample(args):
    
    # Initialize list to hold target records in exact order
    targets = []

    # Initialize dictionary mapping (chrom, pos) tuples to target indexes
    target_keys = {}

    # Open prepared targets file for reading
    with open(args.targets, "r") as tf:
        
        # Iterate through each canonical target record with line index
        for idx, line in enumerate(tf):
            
            # Parse tab-delimited target fields
            c, p, m, r, a = line.rstrip("\n").split("\t")
            
            # Parse position as integer
            p = int(p)
            
            # Append target record tuple to targets list
            targets.append((c, p, m, r, a))
            
            # Define lookup key tuple
            k = (c, p)
            
            # Initialize candidate list for coordinate if not present
            if k not in target_keys:
                
                # Create empty list for coordinate key
                target_keys[k] = []
            
            # Append target index and allele metadata
            target_keys[k].append((idx, m, r, a))

    # Read raw sample VCF stream from zip entry or gzip file
    if args.source_type == "ZIP":
        
        # Open container zip archive in read mode
        with zipfile.ZipFile(args.container, "r") as z:
            
            # Open inner compressed VCF entry from zip archive
            with z.open(args.inner) as raw, gzip.open(raw, "rt") as gz:
                
                # Read all lines from sample VCF stream
                lines_stream = gz.readlines()

    # Handle standard compressed VCF file input
    else:
        
        # Open compressed VCF file for text reading
        with gzip.open(args.container, "rt") as gz:
            
            # Read all lines from sample VCF stream
            lines_stream = gz.readlines()

    # Initialize dictionary to index sample calls by genomic coordinate
    sample_calls = {}

    # Track default FORMAT subfield count for missing field formatting
    default_fmt_len = 1

    # Iterate through all lines in the sample VCF stream
    for line in lines_stream:
        
        # Skip header lines starting with hash mark
        if line.startswith("#"):
            
            # Continue to next line
            continue
        
        # Split data line into tab-delimited columns
        parts = line.rstrip("\n").split("\t")
        
        # Ensure row contains at least 10 columns
        if len(parts) < 10:
            
            # Continue to next line
            continue
        
        # Extract chromosome name
        c = parts[0]
        
        # Extract genomic coordinate as integer
        p = int(parts[1])
        
        # Extract variant identifier
        vid = parts[2]
        
        # Extract reference allele in uppercase
        ref = parts[3].upper()
        
        # Extract alternate alleles string in uppercase
        alt_s = parts[4].upper()
        
        # Split FORMAT specification string by colons
        fmt = parts[8].split(":")
        
        # Update default FORMAT subfield count
        default_fmt_len = len(fmt)
        
        # Split sample values string by colons
        sample_val = parts[9].split(":")

        # Construct call entry dictionary
        call_entry = {
            "vid": vid,
            "ref": ref,
            "alts": alt_s.split(",") if alt_s != "." else [],
            "fmt": fmt,
            "subfields": sample_val,
            "gt": sample_val[0] if sample_val else "./.",
        }

        # Define genomic coordinate key tuple
        k = (c, p)

        # Initialize list for coordinate key if not present
        if k not in sample_calls:
            
            # Create empty list for coordinate key
            sample_calls[k] = []

        # Append call entry to coordinate candidates
        sample_calls[k].append(call_entry)

    # Define missing call representation depending on format mode
    missing_call = "./." + ":." * (default_fmt_len - 1) if args.full_format else "./."

    # Initialize list to hold output column values for each target
    col_values = []

    # Iterate through each canonical target locus in order
    for t_idx, (c, p, mid, t_ref, t_alt) in enumerate(targets):
        
        # Initialize target match variable
        target_match = None

        # Check if sample contains any variant calls at this coordinate
        if (c, p) in sample_calls:
            
            # Retrieve candidate calls at this coordinate
            candidates = sample_calls[(c, p)]

            # Iterate through candidate calls
            for cand in candidates:
                
                # Extract candidate reference allele
                cand_ref = cand["ref"]
                
                # Extract candidate alternate alleles list
                alts = cand["alts"]
                
                # Extract candidate raw genotype string
                raw_gt = cand["gt"]

                # Check if call is homozygous reference
                if raw_gt == "0/0":
                    
                    # Verify reference allele matches target
                    if cand_ref == t_ref or trim_minimal(p, cand_ref, t_ref)[1] == t_ref:
                        
                        # Set target match to homozygous reference
                        target_match = ("0/0", cand)
                        
                        # Break candidate loop
                        break

                # Check if call is uncalled or missing
                elif raw_gt in ("./.", "."):
                    
                    # Set target match to missing
                    target_match = ("./.", cand)
                    
                    # Break candidate loop
                    break

                # Handle alternate allele calls
                else:
                    
                    # Iterate through candidate alternate alleles with index
                    for alt_idx, alt in enumerate(alts, 1):
                        
                        # Check for exact string match with target alleles
                        if cand_ref == t_ref and alt == t_alt:
                            
                            # Set target match to designated alternate index
                            target_match = (f"ALT_{alt_idx}", cand)
                            
                            # Break alternate loop
                            break
                        
                        # Apply minimal parsimony trimming to candidate alleles
                        tp, tr, ta = trim_minimal(p, cand_ref, alt)
                        
                        # Check if trimmed alleles match target alleles exactly
                        if tp == p and tr == t_ref and ta == t_alt:
                            
                            # Set target match to designated alternate index
                            target_match = (f"ALT_{alt_idx}", cand)
                            
                            # Break alternate loop
                            break

                    # Break candidate loop if a match was found
                    if target_match:
                        
                        # Break candidate loop
                        break

        # Handle loci with no matching candidate call
        if target_match is None:
            
            # Append missing call representation
            col_values.append(missing_call)

        # Handle matched candidate call
        else:
            
            # Unpack matched call type and candidate dictionary
            call_type, cand = target_match
            
            # Extract raw genotype string
            raw_gt = cand["gt"]
            
            # Extract sample subfield values list
            sample_val = cand["subfields"]

            # Determine genotype delimiter character (slash or pipe)
            sep = "/" if "/" in raw_gt else "|"
            
            # Split raw genotype into individual allele tokens
            gt_alleles = raw_gt.split(sep)

            # Handle homozygous reference output
            if call_type == "0/0":
                
                # Set output genotype to 0/0
                gt_out = f"0{sep}0"

            # Handle missing call output
            elif call_type == "./.":
                
                # Set output genotype to ./.
                gt_out = f".{sep}."

            # Handle alternate allele call output
            else:
                
                # Extract target alternate allele index string
                target_alt_idx = str(call_type.split("_")[1])
                
                # Initialize validity flag
                is_valid = True
                
                # Initialize list for projected allele tokens
                new_alleles = []

                # Iterate through raw genotype allele tokens
                for a in gt_alleles:
                    
                    # Map reference allele to 0
                    if a == "0":
                        
                        # Append reference allele token
                        new_alleles.append("0")
                    
                    # Map matched target alternate allele to 1
                    elif a == target_alt_idx:
                        
                        # Append alternate allele token
                        new_alleles.append("1")
                    
                    # Flag off-target alleles as invalid
                    else:
                        
                        # Set validity flag to false
                        is_valid = False
                        
                        # Break allele loop
                        break

                # Assign projected genotype string or mask to missing
                gt_out = sep.join(new_alleles) if is_valid else f".{sep}."

            # Handle GT-only output mode
            if not args.full_format:
                
                # Append genotype string to column values
                col_values.append(gt_out)

            # Handle full TVC FORMAT output mode
            else:
                
                # Check if projected call resulted in missing genotype
                if gt_out == f".{sep}.":
                    
                    # Append missing call representation
                    col_values.append(missing_call)
                
                # Emit full format fields with updated genotype
                else:
                    
                    # Create copy of sample subfields list
                    sub_copy = list(sample_val)
                    
                    # Replace first subfield with projected genotype
                    sub_copy[0] = gt_out
                    
                    # Append colon-joined subfield string to column values
                    col_values.append(":".join(sub_copy))

    # Open output column file for writing
    with open(args.output, "w") as out_f:
        
        # Write designated line name as the first header row
        out_f.write(args.line_name + "\n")
        
        # Iterate through all projected target values
        for val in col_values:
            
            # Write target value row
            out_f.write(val + "\n")


# ==============================================================================
# Subcommand 4: assemble-matrix
# ==============================================================================
def cmd_assemble_matrix(args):
    
    # Initialize list to hold open sample file handles
    sample_files = []

    # Initialize list to hold sample column header names
    sample_names = []

    # Iterate through all sample column files by numeric index
    for i in range(1, args.total_samples + 1):
        
        # Define path to sample column TSV file
        cpath = os.path.join(args.cols_dir, f"col_{i}.tsv")
        
        # Open sample column file for reading
        f = open(cpath, "r")
        
        # Read designated sample name from line 1
        sname = f.readline().rstrip("\r\n")
        
        # Append sample name to header list
        sample_names.append(sname)
        
        # Append file handle to active file list
        sample_files.append(f)

    # Define FORMAT column descriptor string based on format mode
    fmt_string = "GT:GQ:DP:FDP:RO:FRO:AO:FAO:SAR:SAF:SRF:SRR" if args.full_format else "GT"

    # Open canonical targets file and output matrix file
    with open(args.targets, "r") as tf, open(args.output, "w") as out_f:
        
        # Construct standard VCF 4.2 data header column list
        header_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + sample_names
        
        # Write tab-delimited #CHROM header line
        out_f.write("\t".join(header_cols) + "\n")

        # Iterate through canonical target definitions
        for tline in tf:
            
            # Parse tab-delimited target locus fields
            c, p, m, r, a = tline.rstrip("\r\n").split("\t")
            
            # Initialize VCF row with fixed standard metadata columns
            row = [c, p, m, r, a, ".", "PASS", ".", fmt_string]
            
            # Read projected call from each open sample column file
            for sf in sample_files:
                
                # Read next call line from sample file
                val = sf.readline().rstrip("\r\n")
                
                # Append sample call or default missing value
                row.append(val if val else "./.")
            
            # Write complete horizontal VCF row
            out_f.write("\t".join(row) + "\n")

    # Close all sample column file handles
    for sf in sample_files:
        
        # Close file handle
        sf.close()


# ==============================================================================
# Subcommand 5: sort-vcf
# ==============================================================================
def cmd_sort_vcf(args):
    
    # Load custom contig ordering from file if provided by user
    custom_contigs = load_contig_order(args.contigs_file) if getattr(args, "contigs_file", None) else []

    # Open input VCF and output sorted VCF files
    with open(args.input, "r") as f_in, open(args.output, "w") as f_out:
        
        # Initialize list to hold data rows for sorting
        data_lines = []

        # Initialize set to collect contig names appearing in input VCF
        seen_contigs = set()

        # Iterate through all lines in the input VCF
        for line in f_in:
            
            # Write header lines directly to the output VCF
            if line.startswith("#"):
                
                # Write header line
                f_out.write(line)
            
            # Collect data lines for coordinate sorting
            else:
                
                # Split data line into tab-delimited fields
                parts = line.rstrip("\n").split("\t")
                
                # Append row parts list to data lines
                data_lines.append(parts)
                
                # Record contig name in seen contigs set
                seen_contigs.add(parts[0])

        # Build contig ranking dictionary based on custom file or seen contigs
        contig_rank, _ = build_contig_rank(custom_contigs, seen_contigs)

        # Sort data rows by contig rank, natural sorting, and coordinate position
        data_lines.sort(key=lambda p: sort_key_record(p[0], int(p[1]), contig_rank))

        # Write sorted data rows to the output VCF
        for parts in data_lines:
            
            # Write tab-delimited data row
            f_out.write("\t".join(parts) + "\n")


# ==============================================================================
# Subcommand 6: generate-qc
# ==============================================================================
def cmd_generate_qc(args):
    
    # Initialize total homozygous reference call counter
    tot_hom_ref = 0

    # Initialize total heterozygous call counter
    tot_het = 0

    # Initialize total homozygous alternate call counter
    tot_hom_alt = 0

    # Initialize total missing call counter
    tot_missing = 0

    # Initialize total evaluated target sites counter
    tot_sites = 0

    # Initialize sample count variable
    n_samples = 0

    # Open sorted VCF and companion QC TSV files
    with open(args.vcf, "r") as vf, open(args.qc_output, "w") as qf:
        
        # Write companion QC TSV header line
        qf.write(
            "#CHROM\tPOS\tID\tREF\tALT\tTOTAL_SAMPLES\tHOM_REF_COUNT\tHET_COUNT\tHOM_ALT_COUNT\tMISSING_COUNT\tCALL_RATE\n"
        )

        # Iterate through all lines in the sorted VCF
        for line in vf:
            
            # Skip header lines
            if line.startswith("#"):
                
                # Continue to next line
                continue
            
            # Parse tab-delimited row columns
            parts = line.rstrip("\n").split("\t")
            
            # Extract locus coordinate and allele fields
            chrom, pos, mid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            
            # Extract sample call columns starting from column index 9
            samples = parts[9:]
            
            # Record number of multiplexed samples
            n_samples = len(samples)

            # Initialize per-marker homozygous reference counter
            hr = 0

            # Initialize per-marker heterozygous counter
            het = 0

            # Initialize per-marker homozygous alternate counter
            ha = 0

            # Initialize per-marker missing counter
            missing = 0

            # Iterate through each sample genotype call
            for s in samples:
                
                # Extract GT subfield
                gt = s.split(":")[0]

                # Tally homozygous reference calls
                if gt in ("0/0", "0|0"):
                    
                    # Increment homozygous reference counter
                    hr += 1

                # Tally heterozygous calls
                elif gt in ("0/1", "1/0", "0|1", "1|0"):
                    
                    # Increment heterozygous counter
                    het += 1

                # Tally homozygous alternate calls
                elif gt in ("1/1", "1|1"):
                    
                    # Increment homozygous alternate counter
                    ha += 1

                # Tally missing calls
                else:
                    
                    # Increment missing counter
                    missing += 1

            # Accumulate global call counts
            tot_hom_ref += hr

            # Accumulate global heterozygous counts
            tot_het += het

            # Accumulate global homozygous alternate counts
            tot_hom_alt += ha

            # Accumulate global missing counts
            tot_missing += missing

            # Increment total sites counter
            tot_sites += 1

            # Compute number of non-missing calls at this site
            called = hr + het + ha

            # Compute marker call rate string formatted to 4 decimals
            call_rate = f"{(called / n_samples):.4f}" if n_samples > 0 else "0.0000"

            # Write marker audit row to companion QC report
            qf.write(f"{chrom}\t{pos}\t{mid}\t{ref}\t{alt}\t{n_samples}\t{hr}\t{het}\t{ha}\t{missing}\t{call_rate}\n")

    # Compute overall total calls across all sites and samples
    total_calls = tot_hom_ref + tot_het + tot_hom_alt + tot_missing

    # Print blank line before console audit summary
    print("")

    # Print top border of console audit summary
    print("================================================================================")

    # Print pipeline completion banner title
    print(" AgriSeq Quick Pipeline Complete!")

    # Print divider line of console audit summary
    print("================================================================================")

    # Print final output VCF path
    print(f" Output VCF File : {args.final_vcf}")

    # Print CSI index note explaining Tabix omission
    print(f" CSI Index File  : {args.final_vcf}.csi (Tabix .tbi omitted due to >512 Mbp wheat chromosomes)")

    # Print companion QC TSV path
    print(f" Companion QC TSV: {args.qc_output}")

    # Print total evaluated sites count
    print(f" Total Sites     : {tot_sites}")

    # Print total multiplexed samples count
    print(f" Total Samples   : {n_samples}")

    # Print formatted genotype call counts
    print(
        f" Genotype Counts : Hom-Ref(0/0)={tot_hom_ref:,} | Het(0/1)={tot_het:,} | "
        f"Hom-Alt(1/1)={tot_hom_alt:,} | Missing(./.)={tot_missing:,}"
    )

    # Check if total calls is greater than zero
    if total_calls > 0:
        
        # Calculate overall percentage call rate
        overall_call_rate = ((total_calls - tot_missing) / total_calls) * 100
        
        # Print overall percentage call rate
        print(f" Overall Call Rate: {overall_call_rate:.2f}%")

    # Print bottom border of console audit summary
    print("================================================================================")


# ==============================================================================
# CLI Dispatcher
# ==============================================================================
def main():
    
    # Initialize argument parser with pipeline description
    parser = argparse.ArgumentParser(description="AgriSeq Quick Multiplexer Helper")

    # Initialize subparser container for subcommands
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: index-samples
    p_idx = subparsers.add_parser("index-samples")
    
    # Add argument for index manifest output path
    p_idx.add_argument("--output", required=True)
    
    # Add argument for input zip directory
    p_idx.add_argument("--zip-dir", default=None)
    
    # Add argument for input pre-extracted directory
    p_idx.add_argument("--input-dir", default=None)
    
    # Add argument for sample mapping file
    p_idx.add_argument("--map-file", default=None)
    
    # Add argument for lines subset list file
    p_idx.add_argument("--lines-file", default=None)
    
    # Set default execution function for index-samples
    p_idx.set_defaults(func=cmd_index_samples)

    # Subcommand: prepare-manifest
    p_man = subparsers.add_parser("prepare-manifest")
    
    # Add argument for raw target manifest input path
    p_man.add_argument("--manifest", required=True)
    
    # Add argument for clean sorted targets output path
    p_man.add_argument("--output-targets", required=True)
    
    # Add argument for ordered contigs output path
    p_man.add_argument("--output-contigs", required=True)
    
    # Add argument for optional headless custom contigs sort order file
    p_man.add_argument("--contigs-file", default=None)
    
    # Set default execution function for prepare-manifest
    p_man.set_defaults(func=cmd_prepare_manifest)

    # Subcommand: project-sample
    p_proj = subparsers.add_parser("project-sample")
    
    # Add argument for source input type (ZIP or FILE)
    p_proj.add_argument("--source-type", required=True, choices=["ZIP", "FILE"])
    
    # Add argument for container file or archive path
    p_proj.add_argument("--container", required=True)
    
    # Add argument for inner entry path inside container zip
    p_proj.add_argument("--inner", default="")
    
    # Add argument for clean targets manifest path
    p_proj.add_argument("--targets", required=True)
    
    # Add argument for single-column output file path
    p_proj.add_argument("--output", required=True)
    
    # Add argument for designated sample line name
    p_proj.add_argument("--line-name", required=True)
    
    # Add flag to retain full TVC FORMAT fields
    p_proj.add_argument("--full-format", action="store_true")
    
    # Set default execution function for project-sample
    p_proj.set_defaults(func=cmd_project_sample)

    # Subcommand: assemble-matrix
    p_asm = subparsers.add_parser("assemble-matrix")
    
    # Add argument for clean targets manifest path
    p_asm.add_argument("--targets", required=True)
    
    # Add argument for directory containing sample column files
    p_asm.add_argument("--cols-dir", required=True)
    
    # Add argument for total sample count integer
    p_asm.add_argument("--total-samples", type=int, required=True)
    
    # Add argument for horizontal matrix output path
    p_asm.add_argument("--output", required=True)
    
    # Add flag to retain full TVC FORMAT fields
    p_asm.add_argument("--full-format", action="store_true")
    
    # Set default execution function for assemble-matrix
    p_asm.set_defaults(func=cmd_assemble_matrix)

    # Subcommand: sort-vcf
    p_sort = subparsers.add_parser("sort-vcf")
    
    # Add argument for unsorted VCF input path
    p_sort.add_argument("--input", required=True)
    
    # Add argument for sorted VCF output path
    p_sort.add_argument("--output", required=True)
    
    # Add argument for optional headless custom contigs sort order file
    p_sort.add_argument("--contigs-file", default=None)
    
    # Set default execution function for sort-vcf
    p_sort.set_defaults(func=cmd_sort_vcf)

    # Subcommand: generate-qc
    p_qc = subparsers.add_parser("generate-qc")
    
    # Add argument for sorted VCF input path
    p_qc.add_argument("--vcf", required=True)
    
    # Add argument for companion QC TSV output path
    p_qc.add_argument("--qc-output", required=True)
    
    # Add argument for final compressed VCF path
    p_qc.add_argument("--final-vcf", required=True)
    
    # Set default execution function for generate-qc
    p_qc.set_defaults(func=cmd_generate_qc)

    # Parse command line arguments from sys.argv
    args = parser.parse_args()

    # Dispatch to the configured subcommand handler function
    args.func(args)


# Check if script is executed as primary entry point
if __name__ == "__main__":
    
    # Invoke main CLI dispatcher
    main()
