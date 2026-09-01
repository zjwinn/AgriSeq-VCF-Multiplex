# AgriSeq Variant Merging, Subsetting, and Genomic Sorting Pipeline

A suite of shell and Python scripts for standardizing, merging, subsetting, and sorting targeted genotyping Variant Call Format (VCF) files produced by AgriSeq sequencing platforms.

---

## Key Features

### 1. Merging & Standardization (`AgriSeq_Merger.sh`)
- **Marker ID Standardization:** Assigns uniform IDs to unnamed variant records (`ID = .` and `INFO/OID = .`) using a user-specified reference prefix, chromosome contig, and zero-padded 9-digit physical positions (e.g., `CHS21_1A004783927`). Supports custom prefixes for non-native/alien contigs.
- **Hierarchical Batch Merging:** Merges hundreds or thousands of individual single-sample VCFs in blocks (default: 100 files) using `bcftools merge` before combining intermediate outputs into the master dataset, avoiding OS open-file limits.
- **Multiallelic Variant Handling:** Accurately merges variants at identical physical coordinates into multiallelic VCF records.
- **Sample Reheadering:** Intersects sequencing sample IDs with a two-column lookup file (`example_sample_map.txt`) to replace raw instrument barcodes with line/cultivar designations using `bcftools reheader`.
- **Marker ID Deduplication:** Optional `-d` flag removes concatenated redundant marker IDs from merged loci, retaining the primary identifier.

### 2. Non-Destructive Subsetting & Sorting (`AgriSeq_Subset_and_Sort.sh`)
- **Target Coordinate Subsetting (`-t, --targets`):** Subsets specific target loci using a 3-column manifest (`#CHROM`, `POS`, `ID`) or 1-column marker ID list (e.g., `example_target_manifest.txt`), eliminating off-target amplicon artifacts.
- **Line Subsetting (`-l, --lines`):** Subsets specific lines/cultivars using a single-column text file (e.g., `example_lines_list.txt`).
- **Complete Metadata & Allele Preservation:** Preserves 100% of reference genome `REF` alleles, all true `ALT` alleles, contig designations, and deep genotype metrics.
- **TASSEL GUI Compatibility:** Sorts records in strict canonical chromosome order (`Chr1A`..`Chr7D` first, then non-native/alien scaffolds alphabetically) and physical coordinate order, resolving "Position out of order" import exceptions in the TASSEL GUI.
- **Dual Indexing:** Automatically generates both CSI (`.csi`) and Tabix (`.tbi`) indexes.

---

## Repository Structure

```text
AgriSeq-VCF-Multiplex/
├── AgriSeq_Merger.sh                       # Core merging and standardization pipeline
├── AgriSeq_Subset_and_Sort.sh              # Core manifest subsetting and sorting pipeline
├── rename_vcf_markers.py                   # Python helper for variant standardization
├── agriseq_vcf_env.yml                     # Conda environment definition
├── example_sample_map.txt                  # Example 2-column sample reheadering map
├── example_target_manifest.txt             # Example 3-column target manifest (#CHROM POS ID)
├── example_lines_list.txt                  # Example 1-column lines/taxa subset list
├── Non_HPC_Submission_Samples/             # Non-HPC execution wrappers
│   ├── run_example.sh                      # Local wrapper for variant merging
│   └── run_subset_example.sh               # Local wrapper for subsetting/sorting
├── SCINet_Submission_Samples/              # HPC SLURM submission scripts
│   ├── example_ceres_submission.slurm      # SLURM script for merging on Ceres
│   └── example_ceres_subset_submission.slurm # SLURM script for subsetting on Ceres
├── README.md                               # Pipeline documentation
├── LICENSE                                 # License file
└── .gitignore                              # Git ignore rules
```

---

## Environment Setup

The pipeline requires `bcftools` (>=1.15), `htslib` (>=1.15), and `python` (>=3.8).

### 1. Build and Activate Conda Environment

```bash
# Clone repository
git clone https://github.com/zjwinn/AgriSeq-VCF-Multiplex.git
cd AgriSeq-VCF-Multiplex

# Create environment from specification
conda env create -f agriseq_vcf_env.yml

# Activate environment
conda activate agriseq_vcf_env
```

### 2. Manual Environment Creation (Alternative)

```bash
conda create -n agriseq_vcf_env -c bioconda -c conda-forge bcftools htslib python -y
conda activate agriseq_vcf_env
```

---

## Example File Formats

### 1. Target Manifest (`example_target_manifest.txt`)
A tab-delimited file with three columns: `#CHROM`, `POS`, and `ID`.

```text
#CHROM	POS	ID
Chr1A	1169383	CHS21_1A001169383
Chr1A	4783927	CHS21_1A004783927
Chr1B	1237446	CHS21_1B001237446
Chr1D	45625	CHS21_1D000045625
TEL10_7E739933655	255	TEL10_7E739933655
```

### 2. Sample Reheadering Map (`example_sample_map.txt`)
A two-column tab-delimited file mapping raw sequencing sample names to line/cultivar designations.

```text
Sample_A01_Plate1	LINE_001
Sample_B01_Plate1	LINE_002
Sample_C01_Plate1	LINE_003
```

### 3. Lines List (`example_lines_list.txt`)
A single-column text file of line/cultivar names to subset.

```text
LINE_001
LINE_002
LINE_003
```

---

## Command-Line Usage

### Stage 1: Merging & Demultiplexing (`AgriSeq_Merger.sh`)

```bash
./AgriSeq_Merger.sh \
  --zip-dir "raw_zip_files" \
  --ref-index "CHS21_OT" \
  --contigs "Chr1A,Chr1B,Chr1D,Chr2A,Chr2B,Chr2D,Chr3A,Chr3B,Chr3D,Chr4A,Chr4B,Chr4D,Chr5A,Chr5B,Chr5D,Chr6A,Chr6B,Chr6D,Chr7A,Chr7B,Chr7D" \
  --output "Merged_Dataset.vcf.gz" \
  --batch-size 100 \
  --map-file "example_sample_map.txt" \
  --non-native "OT" \
  -d
```

#### Options:
- `-z, --zip-dir` *(Required)*: Directory containing individual `.zip` archives delivered by the sequencing provider.
- `-r, --ref-index` *(Required)*: Reference prefix to assign to standardized marker IDs (e.g., `CHS21_OT`).
- `-c, --contigs` *(Required)*: Comma-separated list of standard reference contigs.
- `-o, --output` *(Required)*: Output file path for the master merged VCF (`.vcf.gz`).
- `-b, --batch-size` *(Optional, default: 100)*: Number of VCFs to merge simultaneously per batch.
- `-m, --map-file` *(Optional)*: Two-column sample ID to line name lookup file for reheadering.
- `-p, --non-native` *(Optional)*: Prefix string to insert for non-reference contigs.
- `-d, --dedup` *(Optional)*: Retains only the primary marker ID at multiallelic sites instead of concatenated strings.

---

### Stage 2: Subsetting & Canonical Sorting (`AgriSeq_Subset_and_Sort.sh`)

```bash
./AgriSeq_Subset_and_Sort.sh \
  --input "Merged_Dataset.vcf.gz" \
  --output "Subset_Filtered_Sorted.vcf.gz" \
  --targets "example_target_manifest.txt" \
  --lines "example_lines_list.txt" \
  --memory "16g"
```

#### Options:
- `-i, --input` *(Required)*: Path to input VCF file.
- `-o, --output` *(Required)*: Path for final sorted output VCF (`.vcf.gz`).
- `-t, --targets` *(Optional)*: 3-column target manifest (`#CHROM POS ID`) or 1-column marker ID list.
- `-l, --lines` *(Optional)*: Single-column text file of sample line names to subset.
- `-m, --memory` *(Optional, default: 16g)*: Memory buffer limit for bcftools sorting.

---

## SLURM Execution on SCINet Ceres

Example SLURM job submission scripts are provided in `SCINet_Submission_Samples/`:

```bash
# Submit variant merging job
sbatch SCINet_Submission_Samples/example_ceres_submission.slurm

# Submit subsetting & sorting job
sbatch SCINet_Submission_Samples/example_ceres_subset_submission.slurm
```
