# PONG2 – KIR Genotype Imputation & Model Training

[![R](https://img.shields.io/badge/R-%3E%3D%204.0-blue?logo=r&logoColor=white)](https://www.r-project.org/)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub
release](https://img.shields.io/github/v/release/NormanLabUCD/PONG2)](https://github.com/NormanLabUCD/PONG2/releases/tag/v1.0.0)

**PONG2** is an R package with C++ acceleration (via Rcpp) for
**high-accuracy imputation** and **training** of Killer-cell
Immunoglobulin-like Receptor (**KIR**) genotypes from SNP array data in
the KIR locus (chromosome 19q13.4).

It is optimized for population genetics, immunogenetics, and large-scale
biobank studies requiring reliable KIR allele calls across diverse
ancestries.

**Main CLI commands:** - `impute` – predict KIR alleles from target
PLINK files - `train` – build a new prediction model from reference
genotypes + known KIR calls

------------------------------------------------------------------------

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [CLI Setup & Quick Start](#cli-setup--quick-start)
- [Usage](#usage)
  - [impute command](#impute-command)
  - [train command](#train-command)
- [Improving Imputation Accuracy](#improving-imputation-accuracy)
- [Input & Output Formats](#input--output-formats)
- [Dependencies & External Tools](#dependencies--external-tools)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Citation](#citation)
- [Contact & Support](#contact--support)

------------------------------------------------------------------------

## Overview

PONG2 enables scalable and accurate KIR genotyping by combining:

- Region-specific PLINK2 preprocessing
- Optional local minimac4 pre-imputation for missing variants
- Supervised allele prediction models tailored to the highly polymorphic
  KIR region

It supports both **hg19** and **hg38** assemblies and is particularly
useful for studying immune response variation, HLA–KIR interactions, and
disease associations in diverse populations.

------------------------------------------------------------------------

## Features

- Multi-ancestry pre-trained models (EUR, AMR, AFR, EAS, SAS)
- Automatic handling of hg19 / hg38 coordinate differences
- Configurable SNP missingness threshold
- Built-in local imputation fallback (`--fill-missing`) using minimac4
- Support for external pre-imputation (e.g. Michigan Imputation Server)
- Multi-threading via `--threads`
- Automatic chunked prediction for large biobank datasets (\>2,000
  samples)
- Force-run mode for low SNP match scenarios
- Clean separation of preprocessing and prediction steps

------------------------------------------------------------------------

## Requirements

**R version:** ≥ 4.0

**Required R packages** (loaded at runtime): - `readr` - `tidyverse` -
`parallel`

**System tools** (must be in PATH):

| Tool          | Version | Required                                     |
|---------------|---------|----------------------------------------------|
| PLINK2        | ≥ 2.0   | Always                                       |
| minimac4      | ≥ 4.1.6 | Only with `--fill-missing`                   |
| bgzip & tabix | HTSlib  | Only with `--fill-missing`                   |
| Eagle2        | ≥ 2.4   | Only for pre-phasing before `--fill-missing` |

------------------------------------------------------------------------

## Installation

### From GitHub (recommended)

``` r
# Install remotes if needed
if (!require("remotes", quietly = TRUE)) install.packages("remotes")

# Install PONG2 from GitHub
remotes::install_github("NormanLabUCD/PONG2")
```

### From source tarball

Download
[PONG2_1.0.0.tar.gz](https://github.com/NormanLabUCD/PONG2/releases/download/v1.0.0/PONG2_1.0.0.tar.gz)
then:

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

# Add to PATH (add this to your ~/.bashrc or ~/.bash_profile)
export PATH="$(dirname $PONG2_BIN):$PATH"

# Verify installation
pong2 version
```

------------------------------------------------------------------------

## Usage

``` bash
pong2 <command> [options]
```

### Help

``` bash
pong2 --help              # General overview + list of commands
pong2 --help impute       # Detailed help for imputation
pong2 --help train        # Detailed help for training
pong2 version             # Show version number
```

------------------------------------------------------------------------

### impute command

Predict KIR alleles from a target PLINK dataset.

``` bash
pong2 impute [options]
```

#### Required flags

| Flag             | Description                      | Example              |
|------------------|----------------------------------|----------------------|
| `-i, --bfile`    | PLINK bed/bim/fam prefix (chr19) | `data/chr19`         |
| `-o, --output`   | Output directory                 | `results/imputation` |
| `-l, --locus`    | KIR locus to impute              | `KIR3DL1`            |
| `-a, --assembly` | Genome build                     | `hg19` or `hg38`     |

#### Optional flags

| Flag | Default | Description |
|----|----|----|
| `--filter` | `0.005` | Allele frequency filter threshold (`0.005` or `0.01`) |
| `-t, --threads` | `4` | Number of CPU threads |
| `-f, --force` | `false` | Proceed even if SNP matching rate is low (\<50%) |
| `--fill-missing` | `false` | Impute missing SNPs locally with minimac4 (requires `--vcf`) |
| `--vcf` | — | Pre-phased VCF file required when using `--fill-missing` |

#### Examples

``` bash
# Basic imputation
pong2 impute -i example/chr19 -o output -l KIR3DL1 -a hg19
```

------------------------------------------------------------------------

### train command

Build a new KIR prediction model from reference genotypes and known KIR
calls.

``` bash
pong2 train [options]
```

#### Required flags

| Flag | Description | Example |
|----|----|----|
| `-i, --bfile` | Reference PLINK bed/bim/fam prefix | `data/chr19` |
| `-k, --kfile` | CSV with sample IDs and KIR allele calls | `data/kir_calls.csv` |
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

The KIR file (`--kfile`) must be a CSV with the following structure:

    Sample,KIR3DL1_h1,KIR3DL1_h2
    HG00096,KIR3DL1*001,KIR3DL1*002
    HG00097,KIR3DL1*005,KIR3DL1*015

#### Example

``` bash
pong2 train -i example/chr19 -k -i example/kir_call.csv
```

------------------------------------------------------------------------

## Improving Imputation Accuracy

> **NOTE:** \### KIR Region SNP Overlap between input data and 1KGP
> Overlap rate is computed between your input data and the 1000 Genomes
> Project (1KGP) reference panel in the KIR region (chr19).
>
> | Overlap Rate | Status  | Action                                     |
> |--------------|---------|--------------------------------------------|
> | ≥ 50%        | ✅ Pass | Proceed with PONG2 directly                |
> | \< 50%       | ⚠️ Fail | Run Eagle2 + minimac4 pre-imputation first |

If your SNP matching rate is below 50%, PONG2 provides two strategies:

### Option A: Local pre-imputation with minimac4 (built-in)

Pre-phase your data with Eagle2, then run PONG2 with `--fill-missing`:

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

# Step 2: Run PONG2 with fill-missing
pong2 impute \
  --vcf chr19.phased.vcf.gz \
  -o output \
  -l KIR3DL1 \
  -a hg19 \
  --fill-missing \
  -t 20
```

> ⚠️ **Note:** A pre-phased VCF (`--vcf`) is required with
> `--fill-missing`.

### Option B: External pre-imputation (recommended for highest accuracy)

Pre-impute your chr19 data using a public imputation server before
running PONG2:

- [Michigan Imputation Server](https://imputationserver.sph.umich.edu/)
- [TOPMed Imputation
  Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/)
  (recommended for diverse populations)

------------------------------------------------------------------------

## Input & Output Formats

### Input

| File | Format | Description |
|----|----|----|
| PLINK bfile | `.bed/.bim/.fam` | Genotype data for chr19 |
| KIR file | `.csv` | Sample IDs + phased KIR allele calls (train only) |
| VCF | `.vcf.gz` (bgzipped + tabixed) | Pre-phased VCF (required with `--fill-missing`) |

### Output

| File                  | Description                                      |
|-----------------------|--------------------------------------------------|
| `KIR/<locus>.csv`     | Predicted KIR alleles per sample                 |
| `KIR/<locus>.RData`   | Full prediction object (alleles + probabilities) |
| `<locus>_model.RData` | Trained model object (train only)                |
| `<locus>_test.RData`  | Test genotypes (train only, when `--split < 1`)  |

------------------------------------------------------------------------

## Dependencies & External Tools

| Tool | Purpose | Install |
|----|----|----|
| PLINK2 | Genotype preprocessing | [plink2](https://www.cog-genomics.org/plink/2.0/) |
| Eagle2 | Pre-phasing for imputation | [Eagle](https://alkesgroup.broadinstitute.org/Eagle/) |
| minimac4 | Local SNP imputation | [minimac4](https://github.com/statgen/Minimac4) |
| bgzip/tabix | VCF compression & indexing | [HTSlib](https://www.htslib.org/) |

------------------------------------------------------------------------

## Troubleshooting

| Error | Likely Cause | Fix |
|----|----|----|
| `--vcf is required with --fill-missing` | No VCF provided | Supply pre-phased VCF with `--vcf` |
| `High missing rate` | SNP overlap \< 50% | Run Eagle2 + minimac4, or use `--force` |
| `No model found for locus` | Unsupported locus or wrong filter | Check locus name and `--filter` value |
| `incorrect number of dimensions` | Too few training samples | Verify sample overlap between KIR and PLINK files |
| `plink2 not found` | Not in PATH | Add plink2 to PATH |

------------------------------------------------------------------------

## License

PONG2 is licensed under the **GNU General Public License v3.0**
(GPL-3.0).

You are free to use, modify, and distribute PONG2, provided that
derivative works are distributed under the same license. See
[LICENSE](https://normanlabucd.github.io/PONG2/LICENSE) or the [GNU
GPL-3.0 page](https://www.gnu.org/licenses/gpl-3.0.en.html) for details.

------------------------------------------------------------------------

## Citation

If you use PONG2 in your research, please cite:

> Sadeeq SA, Leaton LA, Kichula KM, Farias TDJ, Font-Porterias N,
> Pollock NR, the Colorado Center for Personalized Medicine, Collora CE,
> Castelli EC, Gignoux CR, Norman PJ. **PONG 2.0: Allele Imputation for
> the Killer Cell Immunoglobulin-Like Receptors.** *Manuscript in
> preparation*, 2026.

------------------------------------------------------------------------

## Contact & Support

- **GitHub Issues** (preferred for bug reports, feature requests,
  questions): <https://github.com/NormanLabUCD/PONG2/issues>

- **Email** (for collaboration or private inquiries):
  <paul.norman@cuanschutz.edu>

- **Lab / Institution:** [Norman
  Lab](https://medschool.cuanschutz.edu/dbmi/lab-pages/norman-lab) \|
  [University of Colorado Anschutz Medical
  Campus](https://www.cuanschutz.edu/) \| [Department of Biomedical
  Informatics](https://news.cuanschutz.edu/dbmi) \| [Immunology &
  Microbiology](https://medschool.cuanschutz.edu/basic-sciences/departments/immunology-and-microbiology)

We aim to respond to issues and emails within 1–3 business days. Thank
you for using PONG2 — happy KIR analysis! 🧬
