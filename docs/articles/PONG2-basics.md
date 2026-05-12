# PONG2 Basics: Installation, Quick Start, and Core Usage

## Overview

PONG2 enables scalable and accurate KIR genotyping by combining:

- Region-specific PLINK2 preprocessing
- C++-accelerated SNP filtering and matching
- Optional local minimac4 pre-imputation for missing variants
- Supervised allele prediction models tailored to the highly polymorphic
  KIR region
- Automatic chunked prediction for large biobank datasets (\>2,000
  samples)

It supports **hg19** and **hg38** assemblies and is particularly useful
for studying immune response variation, HLA–KIR interactions, and
disease associations in diverse populations.

------------------------------------------------------------------------

## Features

- Multi-ancestry pre-trained models (EUR, AMR, AFR, EAS, SAS)
- Fast C++ backend for SNP matching and quality filtering
- Automatic hg19 / hg38 coordinate handling
- Configurable SNP missingness threshold
- Built-in local imputation fallback (`--fill-missing`) using minimac4
- Support for external pre-imputation (e.g. Michigan Imputation Server)
- Multi-threading via `--threads`
- Automatic chunked prediction for large biobank datasets (\>2,000
  samples)
- Force-run mode for low SNP match scenarios
- Built-in colorful help system (`pong2 --help`)

------------------------------------------------------------------------

## Requirements

**R version:** ≥ 4.0

**Required R packages** (loaded at runtime):

- `PONG2` (this package)
- `readr`
- `tidyverse`
- `parallel`

**System tools** (must be in PATH):

| Tool          | Version | Required When                       |
|---------------|---------|-------------------------------------|
| PLINK2        | ≥ 2.0   | Always                              |
| minimac4      | ≥ 4.1.6 | `--fill-missing` only               |
| bgzip & tabix | HTSlib  | `--fill-missing` only               |
| Eagle2        | ≥ 2.4   | Pre-phasing before `--fill-missing` |

------------------------------------------------------------------------

## Installation

### From GitHub (recommended — latest version)

``` r
# Install remotes if needed
if (!require("remotes", quietly = TRUE)) install.packages("remotes")

# Install PONG2
remotes::install_github("NormanLabUCD/PONG2")
```

### From release tarball

Download
[PONG2_1.0.0.tar.gz](https://github.com/NormanLabUCD/PONG2/releases/download/v1.0.0/PONG2_1.0.0.tar.gz)
from the latest release:

``` bash
# Standard install
R CMD INSTALL PONG2_1.0.0.tar.gz

# Custom library path
R CMD INSTALL --library=/your/custom/path PONG2_1.0.0.tar.gz
```

### CLI Setup

After installation, make the `pong2` script executable and add it to
your PATH:

``` bash
# Locate the pong2 script
PONG2_BIN=$(Rscript -e "cat(system.file('scripts', 'pong2', package='PONG2'))")

# Make executable
chmod +x "$PONG2_BIN"

# Add to PATH (add this line to your ~/.bashrc or ~/.bash_profile)
export PATH="$(dirname $PONG2_BIN):$PATH"
```

### Verify installation

``` r
library(PONG2)
packageVersion("PONG2")
```

``` bash
pong2 version
```

------------------------------------------------------------------------

## Quick Start Examples

### 1. Basic imputation

``` bash
pong2 impute \
  -i data/target_chr19 \
  -o results/basic \
  -l KIR3DL1 \
  -a hg38 \
  -t 16
```

### 2. Imputation with missing SNP fill-in

Pre-phase your data first (see [Pre-phasing
section](#pre-phasing-the-kir-region)), then:

``` bash
pong2 impute \
  --vcf data/chr19.phased.vcf.gz \
  -o results/imputed \
  -l KIR3DL1 \
  -a hg38 \
  --fill-missing \
  -t 20
```

> ⚠️ **Note:** `--vcf` (pre-phased VCF) is the **only input** required
> with `--fill-missing`.  
> PLINK files cannot hold phased haplotype data — the pipeline derives
> everything from the VCF.

### 3. Training a new model

``` bash
pong2 train \
  -i data/reference_chr19 \
  -k data/kir_calls.csv \
  -o models/custom \
  -l KIR3DL1 \
  -a hg19 \
  -t 20
```

------------------------------------------------------------------------

## Core Usage Reference

### Help

``` bash
pong2 --help              # General overview + list of commands
pong2 --help impute       # Detailed help for imputation
pong2 --help train        # Detailed help for training
pong2 version             # Show version number
```

------------------------------------------------------------------------

### `impute` command

``` bash
pong2 impute [options]
```

#### Required flags

| Flag | Description | Example |
|----|----|----|
| `-i, --bfile` | PLINK bed/bim/fam prefix (normal imputation) | `data/chr19` |
| `--vcf` | Pre-phased VCF file (required with `--fill-missing`) | `data/chr19.phased.vcf.gz` |
| `-o, --output` | Output directory (created if it doesn’t exist) | `results/imputation` |
| `-l, --locus` | KIR locus to impute | `KIR3DL1` |
| `-a, --assembly` | Genome build | `hg19` or `hg38` |

> **Note:** `-i` and `--vcf` are mutually exclusive: - Normal
> imputation: use `-i` (PLINK bfile) - `--fill-missing`: use `--vcf`
> only (PLINK derived internally from VCF)

#### Optional flags

| Flag | Default | Description |
|----|----|----|
| `--filter` | `0.005` | Allele frequency filter threshold (`0.005` or `0.01`) |
| `-t, --threads` | `4` | Number of CPU threads |
| `-f, --force` | `false` | Proceed even if SNP matching rate is below 50% |
| `--fill-missing` | `false` | Impute missing SNPs locally with minimac4 (requires `--vcf`) |

------------------------------------------------------------------------

### `train` command

``` bash
pong2 train [options]
```

#### Required flags

| Flag | Description | Example |
|----|----|----|
| `-i, --bfile` | Reference PLINK bed/bim/fam prefix | `data/chr19` |
| `-k, --kfile` | CSV with sample IDs and phased KIR allele calls | `data/kir_calls.csv` |
| `-o, --output` | Directory to save trained model | `models/KIR3DL1` |
| `-l, --locus` | KIR locus to train | `KIR3DL1` |
| `-a, --assembly` | Genome build | `hg19` or `hg38` |

#### Optional flags

| Flag | Default | Description |
|----|----|----|
| `-t, --threads` | `4` | Number of CPU threads |
| `--nclassifier` | `100` | Number of ensemble classifiers |
| `--split` | `0.7` | Train/validation split proportion |
| `--kirmaf` | `0.00` | Minimum KIR allele frequency filter |
| `--mac` | `3` | Minimum allele count for SNPs |
| `-r, --region` | Optimized default | Custom KIR region (e.g. `55281035-55295784`) |

#### KIR file format

The KIR file (`--kfile`) must be a comma-separated CSV:

    Sample,KIR3DL1_h1,KIR3DL1_h2
    HG00096,KIR3DL1*001,KIR3DL1*002
    HG00097,KIR3DL1*005,KIR3DL1*015
    HG00099,KIR3DL1*020,KIR3DL1*00302

------------------------------------------------------------------------

## Pre-phasing the KIR Region

Pre-phasing is **required** before using `--fill-missing`. Use Eagle2 to
phase your chr19 data:

### hg19

``` bash
eagle \
  --bfile=chr19 \
  --geneticMapFile=genetic_map_hg19.txt.gz \
  --outPrefix=chr19.phased \
  --chrom=19 \
  --numThreads=20 \
  --bpStart=55000000 \
  --bpEnd=55400000
```

### hg38

``` bash
eagle \
  --bfile=chr19 \
  --geneticMapFile=genetic_map_hg38.txt.gz \
  --outPrefix=chr19.phased \
  --chrom=19 \
  --numThreads=20 \
  --bpStart=54000000 \
  --bpEnd=55000000
```

> Eagle2 outputs a phased VCF (`chr19.phased.vcf.gz`) which is passed
> directly to `--vcf`.

------------------------------------------------------------------------

## Improving Imputation Accuracy

> **NOTE: KIR Region SNP Overlap between input data and 1KGP**  
> Overlap rate is computed between your input data and the 1000 Genomes
> Project (1KGP) reference panel in the KIR region.
>
> | Overlap Rate | Status  | Action                                     |
> |--------------|---------|--------------------------------------------|
> | ≥ 50%        | ✅ Pass | Proceed with PONG2 directly                |
> | \< 50%       | ⚠️ Fail | Run Eagle2 + minimac4 pre-imputation first |

### Option A: Local pre-imputation (built-in, quick)

``` bash
# Step 1: Pre-phase with Eagle2
eagle \
  --bfile=chr19 \
  --geneticMapFile=genetic_map_hg19.txt.gz \
  --outPrefix=chr19.phased \
  --chrom=19 \
  --numThreads=20 \
  --bpStart=55000000 \
  --bpEnd=55400000

# Step 2: Run PONG2 with --fill-missing (VCF only — no -i needed)
pong2 impute \
  --vcf chr19.phased.vcf.gz \
  -o results/imputed \
  -l KIR3DL1 \
  -a hg19 \
  --fill-missing \
  -t 20
```

### Option B: External pre-imputation (recommended for highest accuracy)

Pre-impute your chr19 data using a public server before running PONG2:

**Step 1:** Phase chr19 with Eagle2 (see above)

**Step 2:** Upload phased VCF to [Michigan Imputation
Server](https://imputationserver.sph.umich.edu/) or
[TOPMed](https://imputation.biodatacatalyst.nhlbi.nih.gov/) (recommended
for diverse populations)

- Reference Panel: TOPMed r5
- Chromosome: 19 only

**Step 3:** Download imputed VCF and convert to PLINK:

``` bash
plink2 \
  --vcf imputed.dose.vcf.gz dosage=DS \
  --make-bed \
  --out imputed_chr19
```

**Step 4:** Run PONG2:

``` bash
pong2 impute \
  -i imputed_chr19 \
  -o results/final \
  -l KIR3DL1 \
  -a hg38 \
  --filter 0.005
```

### Option C: Force imputation (not recommended)

Proceed despite low SNP match rate — use only when you understand the
implications:

``` bash
pong2 impute -i chr19 -o results -l KIR3DL1 -a hg19 --force
```

------------------------------------------------------------------------

## Next Steps

- See vignette
  [PONG2-imputation](https://normanlabucd.github.io/PONG2/articles/PONG2-imputation.html)
  for detailed imputation workflow
- See vignette
  [PONG2-training](https://normanlabucd.github.io/PONG2/articles/PONG2-training.html)
  for custom model training
- Report issues: [Open a GitHub
  issue](https://github.com/NormanLabUCD/PONG2/issues/new)
