# AgriSeq Variant Merging, Subsetting, and Genomic Sorting Pipeline

A suite of shell and Python scripts for standardizing, merging, subsetting, and sorting targeted genotyping Variant Call Format (VCF) files produced by AgriSeq sequencing platforms.

---

## Key Features

### 1. Fast Manifest-Driven Multiplexing (`AgriSeq_Merge_Subset_Sort_Quick.sh`) — Recommended
- **Bypasses `bcftools merge` multiallelic & padding issues:** Projects each raw single-sample VCF directly against a canonical 5-column target manifest (`#CHROM POS ID REF ALT`) in parallel. This completely eliminates the multiallelic expansion, coordinate window padding, and header bloat introduced by `bcftools merge`.
- **High-Throughput Parallel Processing:** Processes hundreds or thousands of single-sample VCFs simultaneously using multi-threaded extraction (`-T` / `--threads`), multiplexing an entire 15-plate cohort (1,440 samples) in under 60 seconds.
- **Strict Canonical Biallelic Enforcement:**
  - Preserves genotypes matching canonical manifest alleles (`0/0`, `0/1`, `1/1`).
  - Converts off-target alleles (adjacent amplicon indels, off-target SNPs, STR repeat shifts) and uncalled/missing sites strictly to missing data (`./.`).
- **Flexible FORMAT Support:**
  - Defaults to pure `GT` for ultra-fast, compact TASSEL/T3/GWAS/GS compatibility.
  - Optional `-f` / `--full-format` flag preserves deep TVC FORMAT fields (`GQ`, `DP`, `RO`, `AO`, etc.).
- **Automatic Line Subsetting & Renaming:** Directly filters raw VCFs against an optional lines list (`-l`) and translates raw instrument barcode IDs into cultivar names via a two-column sample map (`-m`).
- **TASSEL GUI Genomic Sorting:** Sorts records in strict canonical chromosome order (`Chr1A`..`Chr7D` first, then non-native/alien scaffolds alphabetically) and physical coordinate order.
- **CSI Indexing (Tabix .tbi Omitted):** Automatically generates a `.csi` index. Tabix (`.tbi`) indexing is intentionally omitted because wheat chromosomes (e.g. `Chr3B` > 512 Mbp) exceed the Tabix index coordinate format limit.
- **Per-Marker QC Audit Report:** Automatically computes and writes a companion summary table (`<output_prefix>_qc.tsv`) reporting per-marker sample counts, call rates, and genotype breakdowns.

### 2. Multi-Step Merging & Standardization (`AgriSeq_Merger.sh`)
- **Marker ID Standardization:** Assigns uniform IDs to unnamed variant records (`ID = .` and `INFO/OID = .`) using a user-specified reference prefix, chromosome contig, and zero-padded 9-digit physical positions (e.g., `CHS21_1A004783927`). Supports custom prefixes for non-native/alien contigs.
- **Hierarchical Batch Merging:** Merges hundreds or thousands of individual single-sample VCFs in blocks (default: 100 files) using `bcftools merge` before combining intermediate outputs into the master dataset, avoiding OS open-file limits.
- **Sample Reheadering:** Intersects sequencing sample IDs with a two-column lookup file (`example_inputs/example_sample_map.txt`) to replace raw instrument barcodes with line/cultivar designations using `bcftools reheader`.
- **Marker ID Deduplication:** Optional `-d` flag removes concatenated redundant marker IDs from merged loci, retaining the primary identifier.

### 3. Subsetting, Biallelic Standardization & Sorting (`AgriSeq_Subset_and_Sort.sh`)
- **Dynamic Target Coordinate + ID Subsetting (`-t, --targets`):** Subsets target loci using a tab-delimited manifest dynamically inspected for required columns (`#CHROM`, `POS`, `ID` by default; `#CHROM`, `POS`, `ID`, `REF`, `ALT` when `-b` is enabled, e.g., `target_manifests/Wheat_AgSeqv2_Targets_Manifest.txt` or `example_inputs/example_target_manifest.txt`). Matching both physical coordinates and marker IDs ensures that only the true designed target variant is retained while rejecting adjacent amplicon indels and off-target (`_OT_`) artifacts.
- **Strict Biallelic / T3 (Breedbase) Compliance Mode (`-b, --biallelic-targets`):** Enforces exact canonical `REF` and `ALT` alleles from the target manifest for T3/Breedbase protocol compliance, GRM construction, GWAS, genomic selection (GS), and MAS:
  - Trims `bcftools merge` 3'-padded windows and 5'-anchored Hotspot blocks back to canonical target coordinates and alleles.
  - Strips off-target secondary alleles (adjacent amplicon indels, off-target SNPs, STR repeat variants) while always preserving the canonical target `ALT` base(s) even at monomorphic `0/0` sites.
  - Dynamically slices `Number=A`, `Number=R`, and `Number=G` `FORMAT`/`INFO` fields to the single target `ALT` allele and masks any sample genotype carrying an off-target allele or half-missing call completely to `./.`.
  - Emits a concise console summary and writes a companion per-marker QC audit table (`<output_prefix>_biallelic_qc.tsv`).
- **Non-Destructive Default Mode:** When `-b` is omitted, preserves 100% of multiallelic `REF`/`ALT` alleles, headers, and deep genotype metrics untouched.
- **Strict Input Sanity Checks:** Dynamically validates manifest headers/columns, verifies positive integer coordinates, flags duplicate coordinates, validates target contigs against VCF headers, verifies sample line overlap, prevents file collisions, and validates VCF integrity. Single-column marker ID lists are explicitly rejected to avoid capturing off-target amplicon indels.
- **Line Subsetting (`-l, --lines`):** Subsets specific lines/cultivars using a single-column text file (e.g., `example_inputs/example_lines_list.txt`).
- **TASSEL GUI Compatibility:** Sorts records in strict canonical chromosome order (`Chr1A`..`Chr7D` first, then non-native/alien scaffolds alphabetically) and physical coordinate order, resolving "Position out of order" import exceptions in the TASSEL GUI.

---

## Repository Structure

```text
AgriSeq-VCF-Multiplex/
├── AgriSeq_Merge_Subset_Sort_Quick.sh        # Fast unified manifest-driven multiplexer (Recommended)
├── AgriSeq_Merger.sh                         # Multi-step hierarchical merging and standardization pipeline
├── AgriSeq_Subset_and_Sort.sh                # Multi-step manifest subsetting, biallelic compliance, and sorting pipeline
├── python_bin/                               # Python helper scripts
│   ├── rename_vcf_markers.py                 # Marker renaming and coordinate padding
│   └── standardize_biallelic_targets.py      # Coordinate+ID verification and strict biallelic REF/ALT standardization
├── conda_yml/                                # Conda environment specifications
│   └── agriseq_vcf_env.yml                   # Conda environment definition
├── target_manifests/                         # Panel target manifests
│   └── Wheat_AgSeqv2_Targets_Manifest.txt    # Canonical 5-column target manifest for Wheat AgriSeq v2 panel
├── example_inputs/                           # Example input templates and lookup tables
│   ├── example_target_manifest.txt           # Example 5-column target manifest (#CHROM POS ID REF ALT)
│   ├── example_sample_map.txt                # Example 2-column sample reheadering map
│   └── example_lines_list.txt                # Example 1-column lines/taxa subset list
├── non_HPC_submission_samples/               # Non-HPC execution wrappers
│   ├── run_quick_example.sh                  # Local wrapper for AgriSeq_Merge_Subset_Sort_Quick.sh
│   ├── run_example.sh                        # Local wrapper for AgriSeq_Merger.sh
│   └── run_subset_example.sh                 # Local wrapper for AgriSeq_Subset_and_Sort.sh
├── SCINet_submission_samples/                # HPC SLURM submission scripts
│   ├── example_ceres_quick_submission.slurm  # SLURM script for AgriSeq_Merge_Subset_Sort_Quick.sh
│   ├── example_ceres_submission.slurm        # SLURM script for AgriSeq_Merger.sh
│   └── example_ceres_subset_submission.slurm # SLURM script for AgriSeq_Subset_and_Sort.sh
├── README.md                                 # Pipeline documentation
├── LICENSE                                   # License file
├── .gitattributes                            # Cross-platform LF line-ending rules
└── .gitignore                                # Git ignore rules
```

---

## Environment Setup

The pipeline requires `bcftools` (>=1.15), `htslib` (>=1.15), and `python` (>=3.8).

Depending on your computing platform, choose either **Environment Modules** (recommended on HPC clusters like Ceres) or **Conda** (recommended for cross-system portability).

### Option A: Environment Modules (Recommended on USDA SCINet Ceres)

Ceres provides pre-compiled, optimized modules for all required tools:

```bash
# Clear any conflicting modules and load required tools
module purge
module load bcftools
module load htslib
module load miniconda3
```

### Option B: Conda Environment (Cross-Platform Stability)

Conda isolates package versions and system libraries, ensuring identical behavior across disparate workstations, servers, and HPC systems:

```bash
# Create environment from specification
conda env create -f conda_yml/agriseq_vcf_env.yml

# Activate environment
conda activate agriseq_vcf_env
```

---

## Example File Formats

### 1. Target Coordinate & Allele Manifest (`target_manifests/Wheat_AgSeqv2_Targets_Manifest.txt` / `example_inputs/example_target_manifest.txt`)
A 5-column tab-delimited file containing `#CHROM`, `POS`, `ID`, `REF`, and `ALT`:

```text
#CHROM	POS	ID	REF	ALT
Chr1A	1169383	CHS21_1A001169383	G	C
Chr1A	4783958	CHS21_1A004783958	G	A
Chr1B	1237446	CHS21_1B001237446	C	G
Chr1D	45625	CHS21_1D000045625	C	T
TEL10_7E739933655	255	TEL10_7E739933655	AG	A
```

### 2. Sample Reheadering Map (`example_inputs/example_sample_map.txt`)
A two-column tab-delimited file mapping raw sequencing sample barcode IDs to line/cultivar designations.

```text
Sample_A01_Plate1	LINE_001
Sample_B01_Plate1	LINE_002
Sample_C01_Plate1	LINE_003
```

### 3. Lines List (`example_inputs/example_lines_list.txt`)
A single-column text file of line/cultivar names to subset.

```text
LINE_001
LINE_002
LINE_003
```

---

## Command-Line Usage

### Fast Unified Pipeline (`AgriSeq_Merge_Subset_Sort_Quick.sh`) — Recommended

Directly multiplex raw single-sample VCFs (from provider zip archives or pre-extracted folders), rename samples, subset lines, enforce canonical manifest alleles, sort canonically for TASSEL, and create CSI indexes in one command:

```bash
./AgriSeq_Merge_Subset_Sort_Quick.sh \
  --zip-dir "raw_zip_files" \
  --targets "target_manifests/Wheat_AgSeqv2_Targets_Manifest.txt" \
  --map-file "example_inputs/example_sample_map.txt" \
  --lines "example_inputs/example_lines_list.txt" \
  --output "Multiplexed_Biallelic.vcf.gz" \
  --threads 16
```

#### Options:
- `-z, --zip-dir`: Directory containing individual `.zip` archives delivered by the sequencing provider.
- `-i, --input-dir`: Directory containing pre-extracted `.vcf.gz` files or plate folders (alternative to `-z`).
- `-t, --targets` *(Required)*: Canonical 5-column target manifest (`#CHROM POS ID REF ALT`).
- `-o, --output` *(Required)*: Final sorted output VCF path (`.vcf.gz`).
- `-m, --map-file` *(Optional)*: Two-column sample ID to line name lookup file for reheadering.
- `-l, --lines` *(Optional)*: Single-column text file of sample line names to subset. If omitted, all discovered samples are included.
- `-f, --full-format` *(Optional)*: Retain deep TVC FORMAT fields (`GQ`, `DP`, `RO`, `AO`, etc.). Defaults to GT-only for speed and compactness.
- `-T, --threads` *(Optional)*: Number of parallel worker threads (default: auto-detected CPU cores).

---

### Alternative Multi-Step Workflow

#### Step 1: Merging & Demultiplexing (`AgriSeq_Merger.sh`)
```bash
./AgriSeq_Merger.sh \
  --zip-dir "raw_zip_files" \
  --ref-index "CHS21_OT" \
  --contigs "Chr1A,Chr1B,Chr1D,Chr2A,Chr2B,Chr2D,Chr3A,Chr3B,Chr3D,Chr4A,Chr4B,Chr4D,Chr5A,Chr5B,Chr5D,Chr6A,Chr6B,Chr6D,Chr7A,Chr7B,Chr7D" \
  --output "Merged_Dataset.vcf.gz" \
  --batch-size 100 \
  --map-file "example_inputs/example_sample_map.txt" \
  --non-native "OT" \
  -d
```

#### Step 2: Subsetting, Biallelic Standardization & Sorting (`AgriSeq_Subset_and_Sort.sh`)
```bash
./AgriSeq_Subset_and_Sort.sh \
  --input "Merged_Dataset.vcf.gz" \
  --output "Subset_Filtered_Sorted.vcf.gz" \
  --targets "target_manifests/Wheat_AgSeqv2_Targets_Manifest.txt" \
  --biallelic-targets \
  --lines "example_inputs/example_lines_list.txt" \
  --memory "16g"
```

---

## SLURM Execution on SCINet Ceres

Example SLURM job submission scripts are provided in `SCINet_submission_samples/`:

```bash
# Submit fast unified multiplexing job (Recommended)
sbatch SCINet_submission_samples/example_ceres_quick_submission.slurm

# Submit multi-step variant merging job
sbatch SCINet_submission_samples/example_ceres_submission.slurm

# Submit multi-step subsetting & sorting job
sbatch SCINet_submission_samples/example_ceres_subset_submission.slurm
```
