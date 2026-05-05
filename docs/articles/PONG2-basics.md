# PONG2 Basics: Installation, Quick Start, and Core Usage

## Overview

PONG2 enables scalable and accurate KIR genotyping by combining:

- Region-specific PLINK2 preprocessing  
- C++-accelerated SNP filtering and matching  
- Optional local minimac4 imputation  
- Supervised allele prediction models tailored to the polymorphic KIR
  region

It supports hg19 and hg38 assemblies and is particularly useful for
studying immune response variation, HLA–KIR interactions, and disease
association in diverse populations.

## Features

- Fast C++ backend for SNP matching and quality filtering  
- Automatic hg19 / hg38 coordinate handling  
- Configurable SNP missingness threshold  
- Built-in local imputation fallback (`--fill-missing`)  
- Support for external pre-imputation (Michigan Imputation Server,
  etc.)  
- Multi-threading via `--threads`  
- Force-run mode and missing SNP imputation strategies  
- Built-in colorful help system (`pong2 --help`)

## Requirements

**R version**  
≥ 4.0

**Required R packages** (loaded at runtime):

- PONG2 (this package)
- readr
- tidyverse
- parallel

**System tools** (must be in PATH):

- PLINK 2.0 (required)
- minimac4 ≥ 4.1.6 (required only when using `--fill-missing`)
- bgzip & tabix (HTSlib – usually bundled with minimac4)

## Installation

### From GitHub (recommended – latest version)

``` r
# Install remotes if needed
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install PONG2
remotes::install_github("NormanLabUCD/PONG2")
```

### From release tarball

``` bash
Download PONG2_1.0.0.tar.gz from the latest release
R CMD INSTALL PONG2_1.0.0.tar.gz
#Or to custom library:
R CMD INSTALL --library="path" PONG2_1.0.0.tar.gz
```

### Verify installation

``` r
library(PONG2)
packageVersion("PONG2")
```

## Quick Start Examples

### Basic imputation

``` bash
pong2 impute \
  -i data/target_chr19 \
  -o results/basic \
  -l KIR \
  -a hg38 \
  -t 16
```

### Imputation with missing SNP fallback

``` bash
pong2 impute \
  -i data/vcf_file \
  -o results/imputed \
  -l KIR \
  -a hg38 \
  --fill-missing \
  --threads 32
```

### Training a new model

``` bash
pong2 train \
  -i data/reference_chr19 \
  -k data/kir_calls.csv \
  -o models/custom \
  -l KIR \
  -a hg19 \
  -t 20
```

## Core Usage Reference

### impute command

``` bash
pong2 impute [options]
```

| Flag | Description | Example Value |
|:---|:---|:---|
| `-i, --input` | PLINK bed/bim/fam prefix (should cover KIR locus on chr19) | `data/my_genotypes_chr19` |
| `-o, --output` | Output directory (will be created if it doesn’t exist) | `results/imputation` |
| `-l, --locus` | Target locus (currently only KIR is supported) | `KIR` |
| `-a, --assembly` | Genome build used in the data | `hg19` or `hg38` |

### Optional flags

| Flag | Default | Description |
|:---|:---|:---|
| `--filter` | `0.01` or `0.005` | KIR genotype quality filter threshold used in model training |
| `-t, --threads` | `20` | Number of CPU threads |
| `-f, --force` | `false` | Proceed even if SNP matching rate is low |
| `--fill-missing` | `false` | Impute missing SNPs locally with minimac4 |
| `-m, --model` | `false` | imputation with your own built model |

### train command

``` bash
pong2 train [options]
```

| Flag | Description |
|:---|:---|
| `-i, --bfile` | Reference PLINK bed/bim/fam prefix |
| `-k, --kfile` | CSV with sample IDs and KIR calls (e.g., sample KIR3DL1_h1 KIR3DL1_h2) |
| `-o, --output` | Directory to save trained model |
| `-l, --locus` | KIR |
| `-a, --assembly` | hg19 or hg38 |

### Optional flags

| Flag               | Default    | Description                           |
|:-------------------|:-----------|:--------------------------------------|
| `-t, --threads`    | `20`       | Number of CPU threads                 |
| `--filter`         | `0.01`     | KIR genotype quality filter threshold |
| `p, or --position` | *Optional* | Optional KIR region                   |

## Pre-phasing the KIR region (required step)

Phasing is mandatory for both local and external imputation. \###
Example with Eagle2 (hg19)

``` bash
eagle \
  --bfile=chr19 \
  --geneticMapFile=genetic_map_hg19.txt.gz \
  --outPrefix=phased_chr19 \
  --chrom=19 \
  --numThreads=50 \
  --bpStart=55200000 \
  --bpEnd=55400000
```

### hg38 version

``` bash
eagle \
  --bfile=../chr19 \
  --geneticMapFile=genetic_map_hg38.txt.gz \
  --outPrefix=phased_chr19 \
  --chrom=19 \
  --numThreads=50 \
  --bpStart=54000000 \
  --bpEnd=55000000
```

## Improving Imputation Accuracy

### Local pre-imputation (built-in)

``` bash
pong2 impute ... --fill-missing
```

### External pre-imputation (recommended for best accuracy)

- Phase chr19 (as above)

- Upload phased VCF to Michigan Imputation Server Panel: TOPMed r5

- Chromosome: 19 only

- Download imputed VCF

- Convert to PLINK:

``` bash
plink2 --vcf imputed.dose.vcf.gz dosage=DS --make-bed --out imputed_ready
```

- Run PONG2:

``` bash
pong2 impute -i imputed_ready -o results/final -l KIR -a hg38 --filter 0.005
```

## Next steps

- See vignette
  [“PONG2-imputation”](https://normanlabucd.github.io/PONG2/articles/PONG2-imputation.html)
  for detailed imputation workflow
- See vignette
  [“PONG2-training”](https://normanlabucd.github.io/PONG2/articles/PONG2-training.html)
  for custom model training
- Report issues: [Open a GitHub
  issue](https://github.com/NormanLabUCD/PONG2/issues/new)
