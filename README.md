# AgriSeq Variant Merging Pipeline

The AgriSeq Variant Merging Pipeline standardizes, aggregates, and reheaders targeted genotyping Variant Call Format (VCF) files produced by AgriSeq sequencing platforms.

---

## Features

- **Marker ID Standardization:** Assigns uniform IDs to unnamed variant records (`ID = .` and `INFO/OID = .`) using reference prefix, chromosome, and zero-padded 9-digit physical positions (e.g., `CHS21_1A004783927`). Supports custom prefixes for non-native contigs.
- **Hierarchical Batch Merging:** Merges individual single-sample VCFs in blocks (default: 100 files) using `bcftools merge` before merging intermediate outputs into the final dataset.
- **Multiallelic Variant Handling:** Consolidates variants at identical physical coordinates into multiallelic VCF records via `bcftools merge`.
- **Sample Reheadering:** Intersects sequencing sample IDs with a two-column lookup file to replace raw IDs with line or cultivar designations using `bcftools reheader`.
- **Marker ID Deduplication:** Optional `-d` flag removes concatenated redundant marker IDs from merged loci, retaining the primary identifier.
- **Automated Directory Cleanup:** Uses `mktemp` for temporary workspace allocation and a POSIX `EXIT` trap for directory removal upon completion or error.

---

## Repository Contents

- `AgriSeq_Merger.sh`: Shell script controlling extraction, marker renaming, indexing, batch merging, and reheadering.
- `rename_vcf_markers.py`: Python script for standardized marker ID generation and header modification.
- `agriseq_vcf_env.yml`: Conda environment specification (`bcftools`, `htslib`, `python`).
- `run_example.sh`: Example bash execution wrapper with parameters.
- `example_ceres_submission.slurm`: SLURM batch job script for the USDA ARS SCINet Ceres cluster.
- `.gitignore`: Rules for excluding raw data archives, sample mapping tables, intermediate outputs, and logs.

---

## Environment Configuration

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

To update dependencies after modifying `agriseq_vcf_env.yml`:

```bash
conda env update -f agriseq_vcf_env.yml --prune
```

### 2. Manual Environment Creation (Alternative)

```bash
conda create -n agriseq_vcf_env -c bioconda -c conda-forge bcftools htslib python -y
conda activate agriseq_vcf_env
```

### 3. Verify Executables

```bash
which bcftools bgzip python
bcftools --version
```

---

## SLURM Execution on SCINet Ceres

SLURM compute nodes run non-interactive subshells where Conda initialization functions are not loaded by default. Conda must be sourced from the user installation before environment activation.

### Environment Activation in Batch Scripts

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate agriseq_vcf_env
```

### Example SLURM Script (`example_ceres_submission.slurm`)

```bash
#!/bin/bash
#SBATCH --job-name=AgriSeq_Merger
#SBATCH --output=agriseq_merge_%j.out
#SBATCH --error=agriseq_merge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=medium
#SBATCH --account=breeding_hwwgru
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.name@usda.gov

echo "Job started on $(hostname) at $(date)"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate agriseq_vcf_env

ZIP_DIR="NAG_Zip_Files"
REF_INDEX="CHS21_OT"
CONTIGS="Chr1A,Chr1B,Chr1D,Chr2A,Chr2B,Chr2D,Chr3A,Chr3B,Chr3D,Chr4A,Chr4B,Chr4D,Chr5A,Chr5B,Chr5D,Chr6A,Chr6B,Chr6D,Chr7A,Chr7B,Chr7D"
OUTPUT_VCF="USDA_AgriSeq_2025_2026.vcf.gz"
MAP_FILE="KSM2026-GYT_Sample_Names.txt"
BATCH_SIZE=100
NON_NATIVE_PREFIX="OT"
DEDUP_FLAG="-d"

./AgriSeq_Merger.sh \
  --zip-dir "$ZIP_DIR" \
  --ref-index "$REF_INDEX" \
  --contigs "$CONTIGS" \
  --output "$OUTPUT_VCF" \
  --batch-size "$BATCH_SIZE" \
  --map-file "$MAP_FILE" \
  --non-native "$NON_NATIVE_PREFIX" \
  $DEDUP_FLAG

echo "Job completed at $(date)"
```

### Ceres Submission and Queue Management

```bash
sbatch example_ceres_submission.slurm
squeue -u $USER
```

---

## Command-Line Usage

### Direct Invocation

```bash
./AgriSeq_Merger.sh \
  -z "NAG_Zip_Files" \
  -r "CHS21_OT" \
  -c "Chr1A,Chr1B,Chr1D,Chr2A,Chr2B,Chr2D,Chr3A,Chr3B,Chr3D,Chr4A,Chr4B,Chr4D,Chr5A,Chr5B,Chr5D,Chr6A,Chr6B,Chr6D,Chr7A,Chr7B,Chr7D" \
  -o "USDA_AgriSeq_2025_2026.vcf.gz" \
  -b 100 \
  -m "KSM2026-GYT_Sample_Names.txt" \
  -p "OT" \
  -d
```

### Options

- `-z, --zip-dir` *(Required)*: Directory containing input `.zip` archives of single-sample VCF files.
- `-r, --ref-index` *(Required)*: Reference index string for marker nomenclature prefix (e.g., `CHS21` or `CHS21_OT`).
- `-c, --contigs` *(Required)*: Comma-separated list of target chromosomes/contigs matching the reference genome.
- `-o, --output` *(Required)*: Path for final merged VCF (`.vcf.gz`).
- `-b, --batch-size` *(Optional, default: 100)*: Integer count of VCF files merged per intermediate batch.
- `-m, --map-file` *(Optional)*: Two-column whitespace-delimited file mapping sequencing IDs to line names.
- `-p, --non-native` *(Optional)*: Prefix string for non-reference contigs.
- `-d, --deduplicate-ids` *(Optional)*: Strips concatenated semicolon-delimited variant IDs, retaining the first ID.
- `-h, --help`: Prints usage documentation and exits.

---

## Pipeline Execution Steps

1. **Extraction:** Allocates temporary directory via `mktemp -d` and unpacks `.zip` archives into separate subdirectories.
2. **Standardization:** Runs `rename_vcf_markers.py` on each `.vcf.gz` to convert missing variant IDs (`ID=.`) into coordinate-based names.
3. **Compression & Indexing:** Compresses standardized VCFs with `bgzip` and indexes with `bcftools index -c` (CSI indexing required for chromosome coordinates > 512 Mbp).
4. **Hierarchical Merging:** Segregates indexed files into batches of size `-b`, executes `bcftools merge -l` per batch, then merges intermediate outputs into the final VCF.
5. **Deduplication (Optional):** If `-d` is passed, truncates semicolon-concatenated marker IDs to the first identifier.
6. **Reheadering (Optional):** If `-m` is provided, filters the mapping file against sample IDs in the VCF and reheaders the file using `bcftools reheader`.
7. **Cleanup:** Removes temporary files via shell `EXIT` trap.

---

## Input File Specifications

### Sample Mapping File (`-m, --map-file`)

Two-column whitespace-delimited text file (no header row):

```text
23786    BOB_DOLE
21427    ZENDA
56553    KS_PROVIDENCE
30137    SHOWDOWN
```
