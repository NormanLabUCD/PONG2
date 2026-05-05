# PONG2 Training: Building Custom KIR Prediction Models

## Overview

This vignette provides a complete, step-by-step guide to training custom
KIR prediction models using the train command in PONG2. Training is
useful when:

-You have a new or population-specific reference dataset -You want to
update models with additional samples -You need locus-specific models
(e.g. only KIR3DL1/S1) -You want to experiment with different SNP
quality filters

Prerequisites:

- PLINK 2.0 installed and in PATH
- R ≥ 4.0 with PONG2 installed
- Reference PLINK files (chr19 or KIR region)
- Known KIR allele calls in CSV format

## Input Data Requirements

1.  Reference genotypes (–bfile)

- PLINK bed/bim/fam files containing SNPs in the KIR locus (chr19)
- Recommended: extract region first to reduce size

``` bash
plink2 \
  --bfile full_reference \
  --chr 19 \
  --make-bed \
  --out reference_chr19
```

2.  Known KIR calls (–kfile)

CSV file with sample IDs and allele calls Minimal format with header
optional:

| sample_id   | KIR2DL1_h1 | KIR2DL1_h2 | KIR3DL1_h1 | KIR3DL1_h2 |
|:------------|:-----------|:-----------|:-----------|:-----------|
| `SAMPLE001` | \*00101    | \*00302    | \*00101    | \*00201    |
| `SAMPLE002` | \*00104    | \*0000     | \*00401    | \*00501    |
| `SAMPLE003` | \*009      | \*00601    | \*00302    | \*00101    |

- Columns: two per KIR locus/gene where h1 and h2 represent KIR
  haplotypes
- Alleles: use standard nomenclature (e.g. *001,* 00201)
- null allele: leave as `*0000`

## Run training

``` bash
pong2 train \
  -i reference_chr19 \
  -k kir_calls.csv \
  -o models/custom \
  -l KIR \
  -a hg38 \
  -t 10
```

## After training completes, the output directory (-o) will contain:

- model.rds or model.RData — trained prediction model (main file)
- training_log.txt — summary of training process, metrics, and warnings
- tmp/ — temporary files (can be deleted)

## Using a custom model for imputation

``` bash
pong2 impute \
  -i chr19 \
  -o results/ \
  -l KIR \
  -a hg38 \
  -m models/custom/model.rds
```

*Note:* The `-m, --model` flag is used to specify a custom model file.

## Troubleshooting Training

- “No matching samples” → Ensure sample IDs in –kfile exactly match
  PLINK .fam
- Low SNP coverage → Increase –filter or pre-impute reference panel
- Memory issues → Reduce threads or use a machine with more RAM
- Slow training → Use more –threads or smaller reference dataset
