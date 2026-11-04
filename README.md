# AgriSeq-VCF-Multiplex

A high-performance pipeline for multiplexing, subsetting, and standardizing targeted amplicon genotyping VCF files from AgriSeq platforms.

---

## Rationale & Methodology

Targeted amplicon sequencing panels assay predefined genomic loci with known reference and alternate alleles. Standard multi-sample variant merging utilities (e.g., `bcftools merge`) dynamically combine proximate variants across samples into unified coordinate intervals. When single samples contain flanking insertions, deletions, or off-target substitutions, this merging behavior introduces multiallelic expansion, coordinate window shifts, and multi-nucleotide padding across the entire cohort. Such artifacts disrupt variant nomenclature and lead to incompatibilities in downstream quantitative genetic software (e.g., TASSEL, GAPIT, Breedbase/T3).

`AgriSeq_Merge_Subset_Sort_Quick.sh` resolves this by evaluating raw single-sample VCFs directly against a canonical 5-column target manifest (`#CHROM`, `POS`, `ID`, `REF`, `ALT`) prior to assembly:
1. **Parallel Manifest Projection:** Single-sample calls are indexed and matched against target loci in parallel.
2. **Allele Standardization:** Genotypes matching canonical target alleles are retained (`0/0`, `0/1`, `1/1`). Off-target variants (flanking indels, unlisted SNPs) and uncalled loci are masked to missing data (`./.`).
3. **Streamed Matrix Assembly:** Standardized sample genotype columns are joined into a uniform VCF 4.2 matrix in lockstep with the manifest.
4. **Flexible Contig Sorting:** Genomic records are sorted according to user-specified contig order (`-c`), canonical wheat order (`Chr1A`–`Chr7D`), or natural alphanumeric order, with ascending physical coordinates.
5. **CSI Indexing:** Output VCFs are compressed with `bgzip` and indexed with CSI format exclusively, avoiding coordinate overflow errors on large chromosomes (>512 Mbp) that occur with Tabix (`.tbi`).
6. **Automated Quality Control:** A companion TSV report calculates per-marker call rates and genotype frequencies across the cohort.

---

## Quick Start

```bash
# Clone the repository
git clone https://github.com/zwwinn/AgriSeq-VCF-Multiplex.git
cd AgriSeq-VCF-Multiplex

# Configure Conda environment (for local execution)
conda env create -f conda_yml/agriseq_vcf_env.yml
conda activate agriseq_vcf_env

# Execute pipeline
./AgriSeq_Merge_Subset_Sort_Quick.sh \
  --zip-dir "raw_zip_files" \
  --targets "target_manifests/Wheat_AgSeqv2_Targets_Manifest.txt" \
  --map-file "example_inputs/example_sample_map.txt" \
  --lines "example_inputs/example_lines_list.txt" \
  --contigs "example_inputs/example_contigs_order.txt" \
  --output "Multiplexed_Biallelic.vcf.gz" \
  --threads 16
```

---

## Command-Line Usage

```bash
./AgriSeq_Merge_Subset_Sort_Quick.sh \
  -t <target_manifest> \
  -o <output_vcf.gz> \
  [-z <zip_dir> | -i <input_dir>] \
  [-m <sample_map>] \
  [-l <lines_list>] \
  [-c <contigs_file>] \
  [-f] \
  [-T <threads>]
```

### Options

| Flag | Argument | Description |
| :--- | :--- | :--- |
| `-t, --targets` | `<file>` | **Required.** Canonical 5-column target manifest (`#CHROM POS ID REF ALT`, tab-delimited). |
| `-o, --output` | `<file>` | **Required.** Path for final compressed and indexed VCF (`.vcf.gz`). |
| `-z, --zip-dir` | `<dir>` | Directory containing raw provider `.zip` archives. *(Either `-z` or `-i` required).* |
| `-i, --input-dir` | `<dir>` | Directory containing pre-extracted `.vcf.gz` files. *(Alternative to `-z`).* |
| `-m, --map-file` | `<file>` | Two-column file mapping raw VCF sample barcode IDs to line names (`<barcode_id> <line_name>`). |
| `-l, --lines` | `<file>` | Single-column file of line names to subset. If omitted, all discovered samples are included. |
| `-c, --contigs` | `<file>` | Headless text file specifying contig sort order (one contig per line). Overrides default wheat order. |
| `-f, --full-format` | *none* | Retain TVC FORMAT subfields (`GQ`, `DP`, `RO`, `AO`, etc.). Default is GT-only. |
| `-T, --threads` | `<int>` | Worker threads for parallel execution (default: auto-detected CPU cores). |
| `-h, --help` | *none* | Display help message and exit. |

---

## Input Formats

### 1. Target Manifest
A 5-column tab-delimited table defining canonical target loci and alleles:

```text
#CHROM	POS	ID	REF	ALT
Chr1A	1169383	CHS21_1A001169383	G	C
Chr1A	4783958	CHS21_1A004783958	G	A
Chr1B	1237446	CHS21_1B001237446	C	G
Chr1D	45625	CHS21_1D000045625	C	T
```

### 2. Sample Name Map (Optional)
A 2-column whitespace- or tab-delimited table mapping instrument barcode identifiers to sample/line designations:

```text
23786	BOB_DOLE
21427	ZENDA
56553	KS_PROVIDENCE
```

### 3. Lines List (Optional)
A single-column text file of sample/line names to include. Samples not listed in this file are excluded during initial indexing:

```text
BOB_DOLE
ZENDA
KS_PROVIDENCE
```

### 4. Contig Order File (Optional)
A headless text file specifying the desired chromosome/contig ordering (one per line). Contigs appearing in the manifest but omitted from this file are appended automatically using natural alphanumeric sorting:

```text
Chr1A
Chr1B
Chr1D
Chr2A
...
Chr7D
```

---

## Output Files

1. **Multiplexed VCF (`<output>.vcf.gz`)**: Standard VCF 4.2 file containing biallelic genotypes for all samples across all manifest loci.
2. **CSI Index (`<output>.vcf.gz.csi`)**: Coordinate index for fast querying via `bcftools` and IGV.
3. **Quality Control Report (`<output>_qc.tsv`)**: Tab-delimited per-marker summary:
   - `#CHROM`, `POS`, `ID`, `REF`, `ALT`: Target locus metadata.
   - `TOTAL_SAMPLES`: Cohort sample size.
   - `HOM_REF_COUNT`, `HET_COUNT`, `HOM_ALT_COUNT`: Genotype call counts.
   - `MISSING_COUNT`: Count of uncalled or masked off-target genotypes.
   - `CALL_RATE`: Proportion of non-missing genotypes.

---

## High-Performance Computing (SCINet Ceres)

SLURM-managed HPC environments (e.g., USDA SCINet Ceres) execute compute jobs in isolated shells on compute nodes. Environment configuration (`module purge`, `module load bcftools htslib miniconda3`, and conda activation) is embedded directly within the batch submission script rather than executed on the login node.

To run on Ceres:
1. Edit input paths and SLURM parameters in `SCINet_submission_samples/example_ceres_quick_submission.slurm`.
2. Submit the job from the repository root:

```bash
sbatch SCINet_submission_samples/example_ceres_quick_submission.slurm
```

---

## Repository Structure

```text
AgriSeq-VCF-Multiplex/
├── AgriSeq_Merge_Subset_Sort_Quick.sh        # Primary pipeline execution script
├── python_bin/                               # Modular Python backend
│   └── quick_multiplex.py                    # Indexing, projection, assembly, and QC routines
├── conda_yml/                                # Environment configuration
│   └── agriseq_vcf_env.yml                   # Conda package dependencies
├── target_manifests/                         # Panel target manifests
│   └── Wheat_AgSeqv2_Targets_Manifest.txt    # Canonical 3,999-marker wheat manifest
├── example_inputs/                           # Example input templates
│   ├── example_target_manifest.txt           # Sample manifest format
│   ├── example_sample_map.txt                # Sample barcode-to-name mapping
│   ├── example_lines_list.txt                # Sample filtering list
│   └── example_contigs_order.txt             # Headless contig order file
├── non_HPC_submission_samples/               # Workstation execution scripts
│   └── run_quick_example.sh                  # Local shell wrapper
├── SCINet_submission_samples/                # Cluster submission scripts
│   └── example_ceres_quick_submission.slurm  # SLURM submission template
├── legacy_code/                              # Historical multi-step scripts (pre-manifest projection)
└── README.md                                 # Documentation
```
