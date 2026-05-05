# PONG2 – KIR Genotype Imputation & Model Training

[![R](https://img.shields.io/badge/R-%3E%3D%204.0-blue?logo=r&logoColor=white)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub release](https://img.shields.io/github/v/release/NormanLabUCD/PONG2)](https://github.com/NormanLabUCD/PONG2/releases/tag/v1.0.0)
<!-- Add real CI badge if you set up GitHub Actions later -->
<!-- [![R-CMD-check](https://github.com/NormanLabUCD/PONG2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NormanLabUCD/PONG2/actions) -->

**PONG2** is an R package with C++ acceleration (via Rcpp) for **high-accuracy imputation** and **training** of Killer-cell Immunoglobulin-like Receptor (**KIR**) genotypes from SNP array in the KIR locus (chromosome 19).

It is optimized for population genetics, immunogenetics, and large-scale biobank studies requiring reliable KIR allele calls.

**Main CLI commands:**

- `impute` – predict KIR alleles from target PLINK files  
- `train` – build a new prediction model from reference genotypes + known KIR calls

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
  - [From GitHub (recommended for latest version)](#from-github-recommended-for-latest-version)
  - [From source tarball](#from-source-tarball)
  - [From source directory (developer mode)](#from-source-directory-developer-mode)
- [CLI Setup & Quick Start](#cli-setup--quick-start)
- [Usage](#usage)
  - [impute command](#impute-command)
  - [train command](#train-command)
- [Improving Imputation Accuracy](#improving-imputation-accuracy)
- [Input & Output Formats](#input--output-formats)
- [Dependencies & External Tools](#dependencies--external-tools)
- [Troubleshooting](#troubleshooting)
- [Development & Contributing](#development--contributing)
- [License](#license)
- [Citation](#citation)
- [Contact & Support](#contact--support)

## Overview

PONG2 enables scalable and accurate KIR genotyping by combining:

- region-specific PLINK2 preprocessing  
- optional local minimac4 imputation for missing variants  
- supervised allele prediction models tailored to the polymorphic KIR region

It supports both hg19 and hg38 assemblies and is particularly useful for studying immune response variation, HLA–KIR interactions, and disease association in diverse populations.

## Features
- Automatic handling of hg19 / hg38 coordinate differences  
- Configurable SNP missingness threshold  
- Built-in local imputation fallback (`--fill-missing`)  
- Support for external pre-imputation (e.g. Michigan Imputation Server)  
- Multi-threading via `--threads`  
- Force-run mode and missing SNP imputation strategies  
- Clean separation of preprocessing and prediction steps

## Requirements

**R version**  
≥ 4.0

**Required R packages** (loaded at runtime):
- readr
- tidyverse
- parallel

**System tools** (must be in PATH):

- PLINK 2.0 (required)
- minimac4 ≥ 4.1.6 (required only when using `--fill-missing`)
- bgzip & tabix (HTSlib – usually bundled with minimac4)

## Installation

### From GitHub (recommended for latest version)

Install the development or stable version directly from the repository:

```r
# Install remotes if needed
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install PONG2 from GitHub
remotes::install_github("https://github.com/NormanLabUCD/PONG2")
```

### Standard install
Download the source tarball from the release: [PONG2_1.0.0.tar.gz](https://github.com/NormanLabUCD/PONG2/releases/download/v1.0.0/PONG2_1.0.0.tar.gz)
```bash
R CMD INSTALL PONG2_1.0.0.tar.gz
```

### Custom library path
```bash
R CMD INSTALL --library="installation path" PONG2_1.0.0.tar.gz
```

### USAGE
```bash
pong2 <command> [options]
```
### Basic imputation
```bash
pong2 impute -i bfile -o output -l KIR -a hg19
```

### Imputation with local missing SNP fill-in
```bash
pong2 impute -i bfile -o output -l KIR -a hg19 --fill-missing -t 20
```

### Train a new model
```bash
pong2 train -i bfile -k kfile -o output -l KIR -a hg19 -t 20
```
### Specific detailed help
```bash
pong2 --help                # General overview + list of commands
pong2 --help impute         # Detailed help for imputation
pong2 --help train          # Detailed help for training
pong2 version               # Show version number
```

### impute command
```bash
pong2 impute [options]
```
| Flag | Description | Example Value |
| :--- | :--- | :--- |
| `-i, --input` | PLINK bed/bim/fam prefix (should cover KIR locus on chr19) | `data/my_genotypes_chr19` |
| `-o, --output` | Output directory (will be created if it doesn't exist) | `results/imputation` |
| `-l, --locus` | Target locus (currently only KIR is supported) | `KIR` |
| `-a, --assembly` | Genome build used in the data | `hg19` or `hg38` |

### Optional flags

| Flag | Default | Description |
| :--- | :--- | :--- |
| `--filter` | `0.01` or `0.005` | KIR genotype quality filter threshold used in model training |
| `-t, --threads` | `20` | Number of CPU threads |
| `-f, --force` | `false` | Proceed even if SNP matching rate is low |
| `--fill-missing` | `false` | Impute missing SNPs locally with minimac4 |
| `-m, --model` | `false` | imputation with your own built model |


### train command
```bash
pong2 train [options]
```
| Flag | Description |
| :--- | :--- |
| `-i, --bfile` | Reference PLINK bed/bim/fam prefix |
| `-k, --kfile` | CSV with sample IDs and KIR calls (e.g., sample KIR3DL1_h1 KIR3DL1_h2) |
| `-o, --output` | Directory to save trained model |
| `-l, --locus` | KIR |
| `-a, --assembly` | hg19 or hg38 |

### Optional flags
| Flag | Default | Description |
| :--- | :--- | :--- |
| `-t, --threads` | `20` | Number of CPU threads |
| `--filter` | `0.01` | KIR genotype quality filter threshold |
| `p, or --position` | *Optimized KIR region* | Optional KIR region |

### Examples
#### 1.
```bash
pong2 impute -i example/chr19 -o output -l KIR3DL1 -a hg19
```
#### 2. Pre-imputation (recommended for best accuracy)
If the initial SNP matching rate in the KIR region is low (e.g., < 50%), PONG2 provides two main 
strategies to improve accuracy by imputing missing or low-quality SNPs **before** running the core PONG2 prediction step.

#### Note: Pre-phasing the KIR region is required for both local and external pre-imputation (https://alkesgroup.broadinstitute.org/Eagle/)

```bash
eagle \
  --bfile=chr19 \
  --geneticMapFile=genetic_map_hg19.txt.gz \
  --outPrefix=chr19.phased \
  --chrom=19 \
  --numThreads=50 \
  --bpStart=55200000 \
  --bpEnd=55300000
```
- Option A: Run PONG2 with local pre-imputation using minimac4 (built-in – quick & automated) using flag the `--fill-missing`

```bash
#Use the `--fill-missing` flag when running `impute`:
pong2 impute -i chr19.phased.vcf -o output -l KIR3DL1 -a hg19 --fill-missing -t 32
```

- Option B: External pre-imputation (recommended for highest accuracy)

pre-impute your chr19 data using a public imputation server before running PONG2.
Recommended server: Michigan Imputation Server (https://imputationserver.sph.umich.edu/)

- You may chose to force PONG2 to continue even with low SNP matching rate using the flag `-f or --force`
```bash
pong2 impute -i example/chr19 -o output -l KIR3DL1 -a hg19 -f
```
#### 3.
```bash
pong2 train -i example/chr19 -k example/kir_calls.csv -o output -l KIR -a hg19 -t 24
```
**NOTE:**
### KIR Region SNP Overlap between input data and 1KGP

Overlap rate is computed between your input data and the 1000 Genomes Project (1KGP) 
reference panel in the KIR region.

| Overlap Rate | Status | Action |
|-------------|--------|--------|
| ≥ 50% | ✅ Pass | Proceed with PONG2 |
| < 50% | ⚠️ Fail | Run Eagle2 + minimac4 pre-imputation first |

**Detailed tutorials & examples** are available in the vignettes:  
[PONG2 Basics](https://normanlabucd.github.io/PONG2/) • [Imputation Workflow](https://normanlabucd.github.io/PONG2/articles/PONG2-imputation.html) • [Training Models](https://normanlabucd.github.io/PONG2/articles/PONG2-training.html)


### License



PONG2 is licensed under the **GNU General Public License v3.0** (GPL-3.0).

Full license text: [LICENSE](LICENSE) (included in the repository root)

You are free to use, modify, and distribute PONG2, provided that:

- You distribute any derivative work under the same GPL-3.0 license
- You include the original copyright notice and license text

See the [GNU GPL-3.0 page](https://www.gnu.org/licenses/gpl-3.0.en.html) for details.


### Author 

```markdown
Suraju A Sadeeq, Laura Ann Leaton, Katherine M Kichula, Neus Font-Porterias, Nicholas R Pollock,
the Colorado Center for Personalized Medicine, Christopher E Collora, Erick C. Castelli, Christopher R Gignoux, Paul J Norman
```


### Contact & Support

- **GitHub Issues** (preferred for bug reports, feature requests, questions):  
  [https://github.com/NormanLabUCD/PONG2/issues](https://github.com/NormanLabUCD/PONG2/issues)

- **Email** (for collaboration, private inquiries, or urgent support):  
  paul.norman@cuanschutz.edu

- **Lab / Institution**:  
  [Norman Lab](https://medschool.cuanschutz.edu/dbmi/lab-pages/norman-lab)  
  [University of Colorado Anschutz Medical Campus](https://www.cuanschutz.edu/)  
  [Department of Biomedical Informatics](https://news.cuanschutz.edu/dbmi), [Immunology & Microbiology](https://medschool.cuanschutz.edu/basic-sciences/departments/immunology-and-microbiology)

We aim to respond to issues and emails within 1–3 business days.  
Thank you for using PONG2 — happy KIR analysis! 🧬

