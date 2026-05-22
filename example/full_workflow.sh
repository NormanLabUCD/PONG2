#!/bin/bash
# =================================================================================
# PONG2 Complete Workflow
# From raw SNP array data to KIR genotypes and custom model training and evaluation
# =================================================================================

# =================================================================================
# PART 1: IMPUTATION — Predict KIR alleles from existing models
# =================================================================================

# Step 1: Extract chr19 from full genome
plink2 \
  --bfile full_genome \
  --chr 19 \
  --make-bed \
  --out chr19

# Step 2: Run basic imputation (checks SNP overlap automatically)
pong2 impute \
  -i chr19 \
  -o output/imputation \
  -l KIR3DL1 \
  -a hg19 \
  -t 20

# Step 3a: If SNP overlap < 50% — local pre-imputation with minimac4
# Pre-phase with Eagle2
eagle \
  --bfile=chr19 \
  --geneticMapFile=genetic_map_hg19.txt.gz \
  --outPrefix=chr19.phased \
  --chrom=19 \
  --numThreads=20 \
  --bpStart=55000000 \
  --bpEnd=55400000

# Run PONG2 with local fill-missing (VCF only — no -i needed)
pong2 impute \
  --vcf chr19.phased.vcf.gz \
  -o output/imputation \
  -l KIR3DL1 \
  -a hg19 \
  -t 20 \
  --filter 0.005 \
  --fill-missing

# Step 3b: If SNP overlap < 50% — external pre-imputation (recommended)
# Upload chr19.phased.vcf.gz to Michigan Imputation Server:
#   https://imputationserver.sph.umich.edu/
#   Reference panel: TOPMed r5
#   Chromosome: 19 only
# Download imputed VCF then convert to PLINK
plink2 \
  --vcf imputed.dose.vcf.gz dosage=DS \
  --import-dosage-certainty 0.3 \
  --make-bed \
  --out imputed_chr19

# Run PONG2 on imputed data
pong2 impute \
  -i imputed_chr19 \
  -o output/imputation \
  -l KIR3DL1 \
  -a hg19 \
  -t 20 \
  --filter 0.005

# =============================================================================
# PART 2: TRAINING — Build a custom KIR prediction model
# =============================================================================

# Step 1: Prepare reference genotypes
plink2 \
  --bfile reference_full_genome \
  --chr 19 \
  --make-bed \
  --out reference_chr19

# Step 2: Prepare KIR allele calls CSV
# Format: Sample,KIR3DL1_h1,KIR3DL1_h2
# Example: HG00096,KIR3DL1*001,KIR3DL1*002

# Step 3: Train custom model
pong2 train \
  -i reference_chr19 \
  -k kir_calls.csv \
  -o output/models \
  -l KIR3DL1 \
  -a hg19 \
  -t 20 \
  --nclassifier 100 \
  --split 0.8 \
  --kirmaf 0.005

# Step 4: Evaluate trained model
pong2 evaluate \
  --model-dir output/models \
  --locus KIR3DL1 \
  --threshold 0.5

# Step 5: Use custom model for imputation
pong2 impute \
  -i chr19 \
  -o output/custom_imputation \
  -l KIR3DL1 \
  -a hg19 \
  -t 20 \
  -m output/models/KIR3DL1_model.RData
