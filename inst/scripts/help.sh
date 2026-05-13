#!/usr/bin/env bash
# PONG2 Help System
# KIR Genotype Imputation and Model Training Toolkit
# Version: 1.0.0

set -euo pipefail

# ────────────────────────────────────────────────
# Color & Style Definitions
# ────────────────────────────────────────────────
BOLD=$(tput bold)
GREEN=$(tput setaf 2)
BLUE=$(tput setaf 4)
CYAN=$(tput setaf 6)
YELLOW=$(tput setaf 3)
RED=$(tput setaf 1)
WHITE=$(tput setaf 7)
RESET=$(tput sgr0)

# ────────────────────────────────────────────────
# Main Help Screen
# ────────────────────────────────────────────────
show_main_help() {
    cat << HELP
${BOLD}${GREEN}PONG2 ─ KIR Genotype Imputation & Training Toolkit${RESET}
${BLUE}Version:${RESET}     1.0.0
${BLUE}Description:${RESET} High-accuracy KIR allele imputation and model training from SNP array data (chr19 locus)

${BOLD}${GREEN}USAGE${RESET}
  pong2 <command> [options]

${BOLD}${GREEN}AVAILABLE COMMANDS${RESET}
  ${CYAN}impute${RESET}     Predict KIR alleles from target genotype data
  ${CYAN}train${RESET}      Train a new KIR prediction model from reference data

${BOLD}${GREEN}GLOBAL OPTIONS${RESET}
  ${CYAN}-h, --help [command]${RESET}   Show this help message or command-specific help
  ${CYAN}-v, --version${RESET}          Show version information

${BOLD}${GREEN}GET COMMAND-SPECIFIC HELP${RESET}
  pong2 --help impute
  pong2 --help train

${BOLD}${GREEN}EXAMPLES${RESET}
  # Basic imputation
  pong2 impute -i data/chr19 -o results/ -l KIR3DL1 -a hg38

  # Imputation with local missing SNP fill-in (pre-phased VCF required)
  pong2 impute --vcf chr19.phased.vcf.gz -o results/ -l KIR3DL1 -a hg38 --fill-missing -t 20

  # Force imputation despite low SNP match rate
  pong2 impute -i data/chr19 -o results/ -l KIR3DL1 -a hg19 --force

  # Train a new model
  pong2 train -i data/ref_chr19 -k data/kir_calls.csv -o models/ -l KIR3DL1 -a hg19 -t 20

${BOLD}${YELLOW}IMPORTANT NOTES${RESET}
  • Input must be PLINK format (bed/bim/fam) covering chr19 KIR region
  • For --fill-missing: supply a pre-phased VCF (--vcf) — PLINK cannot hold phased data
  • Pre-phase with Eagle2 before using --fill-missing (see pong2 --help impute)
  • For large datasets (>2,000 samples), prediction is auto-chunked to prevent memory issues
  • Use --fill-missing or external pre-imputation when SNP overlap with 1KGP < 50%

${BOLD}${BLUE}DOCUMENTATION & SUPPORT${RESET}
  GitHub:       https://github.com/NormanLabUCD/PONG2
  Issues:       https://github.com/NormanLabUCD/PONG2/issues
  Vignettes:    https://normanlabucd.github.io/PONG2/
  Citation:     Sadeeq et al. (2026). PONG 2.0: Allele imputation for KIR. Manuscript in preparation.
HELP
}

# ────────────────────────────────────────────────
# Impute-specific detailed help
# ────────────────────────────────────────────────
show_impute_help() {
    cat << HELP
${BOLD}${GREEN}pong2 impute ─ KIR Allele Imputation${RESET}

${BOLD}${WHITE}DESCRIPTION${RESET}
  Predicts KIR genotypes from target SNP array data using pre-trained models.
  Automatically checks SNP overlap with the 1KGP reference panel and routes
  accordingly. For large datasets (>2,000 samples), prediction is automatically
  chunked to prevent memory issues.

${BOLD}${WHITE}REQUIRED OPTIONS${RESET}
  ${CYAN}-i, --bfile FILE${RESET}     PLINK bed/bim/fam prefix — required for normal imputation
  ${CYAN}--vcf FILE${RESET}           Pre-phased VCF file — required when using --fill-missing
  ${CYAN}-o, --output DIR${RESET}     Output directory (created if it does not exist)
  ${CYAN}-l, --locus STR${RESET}      KIR locus to impute (e.g. KIR3DL1, KIR2DL1)
  ${CYAN}-a, --assembly STR${RESET}   Genome build: hg19 or hg38

  ${YELLOW}Note: -i and --vcf are mutually exclusive:${RESET}
    Normal imputation:  use ${CYAN}-i${RESET} (PLINK bfile)
    --fill-missing:     use ${CYAN}--vcf${RESET} only — PLINK derived internally from VCF

${BOLD}${WHITE}OPTIONAL OPTIONS${RESET}
  ${CYAN}-t, --threads INT${RESET}    Number of CPU threads [default: 4]
  ${CYAN}--filter FLOAT${RESET}       Allele frequency filter threshold: 0.005 or 0.01 [default: 0.005]
  ${CYAN}-f, --force${RESET}          Proceed even if SNP match rate is below 50%
  ${CYAN}--fill-missing${RESET}       Impute missing SNPs locally with minimac4 (requires --vcf)

${BOLD}${WHITE}SNP OVERLAP CHECK${RESET}
  PONG2 automatically checks your SNP overlap with the 1KGP reference panel:

  ${GREEN}≥ 50%${RESET}   ✅ Proceed with imputation directly
  ${RED}< 50%${RESET}   ⚠️  Pre-imputation recommended (see below)

${BOLD}${WHITE}PRE-IMPUTATION STRATEGIES${RESET}

  ${BOLD}Option A: Local (built-in) — quick${RESET}
    Step 1: Pre-phase with Eagle2
      eagle --bfile=chr19 --geneticMapFile=genetic_map_hg19.txt.gz \\
            --outPrefix=chr19.phased --chrom=19 --numThreads=20 \\
            --bpStart=55000000 --bpEnd=55400000

    Step 2: Run PONG2 (--vcf only — no -i needed)
      pong2 impute --vcf chr19.phased.vcf.gz -o results/ \\
                   -l KIR3DL1 -a hg19 --fill-missing -t 20

  ${BOLD}Option B: External (recommended for highest accuracy)${RESET}
    Phase chr19 → upload to Michigan Imputation Server → download imputed VCF
    → convert to PLINK → run pong2 impute -i imputed_chr19 ...

  ${BOLD}Option C: Force run (not recommended)${RESET}
    pong2 impute -i chr19 -o results/ -l KIR3DL1 -a hg19 --force

${BOLD}${WHITE}EXAMPLES${RESET}
  # Basic imputation
  pong2 impute -i data/chr19 -o results/ -l KIR3DL1 -a hg38

  # With pre-phased VCF and fill-missing
  pong2 impute --vcf chr19.phased.vcf.gz -o results/ -l KIR3DL1 -a hg38 --fill-missing -t 20

  # Force despite low match rate
  pong2 impute -i data/chr19 -o results/ -l KIR3DL1 -a hg19 --force

${BOLD}${WHITE}OUTPUT FILES${RESET}
  KIR/<locus>.csv          Predicted KIR alleles per sample (main results)
  KIR/<locus>.RData        Full prediction object including allele probabilities

${BOLD}${WHITE}TROUBLESHOOTING${RESET}
  • Low SNP match rate        → Use --fill-missing or pre-impute externally
  • --vcf required error      → Supply pre-phased VCF with --fill-missing
  • minimac4 not found        → Install minimac4 ≥ 4.1.6 and add to PATH
  • 1 sample after --keep     → Check sample ID format (FAM vs KIR file mismatch)
HELP
}

# ────────────────────────────────────────────────
# Train-specific detailed help
# ────────────────────────────────────────────────
show_train_help() {
    cat << HELP
${BOLD}${GREEN}pong2 train ─ Train a New KIR Prediction Model${RESET}

${BOLD}${WHITE}DESCRIPTION${RESET}
  Builds a new KIR prediction model from reference genotypes and known KIR calls.
  Useful for custom populations, updated reference panels, or locus-specific retraining.
  Training is parallelized across CPU threads for efficiency.

${BOLD}${WHITE}REQUIRED OPTIONS${RESET}
  ${CYAN}-i, --bfile FILE${RESET}     PLINK bed/bim/fam prefix (reference genotypes, chr19)
  ${CYAN}-k, --kfile FILE${RESET}     CSV file with sample IDs and phased KIR allele calls
  ${CYAN}-o, --output DIR${RESET}     Directory to save trained model files
  ${CYAN}-l, --locus STR${RESET}      KIR locus to train (e.g. KIR3DL1, KIR2DL1)
  ${CYAN}-a, --assembly STR${RESET}   Genome build: hg19 or hg38

${BOLD}${WHITE}OPTIONAL OPTIONS${RESET}
  ${CYAN}-t, --threads INT${RESET}       Number of CPU threads [default: 4]
  ${CYAN}--nclassifier INT${RESET}       Number of ensemble classifiers [default: 100]
  ${CYAN}--split FLOAT${RESET}           Train/validation split proportion [default: 0.7]
  ${CYAN}--kirmaf FLOAT${RESET}          Minimum KIR allele frequency filter [default: 0.00]
  ${CYAN}--mac INT${RESET}               Minimum allele count for SNPs [default: 3]
  ${CYAN}-r, --region STR${RESET}        Custom KIR SNP region e.g. 55281035-55295784 [default: optimized]

${BOLD}${WHITE}KIR FILE FORMAT${RESET}
  CSV with header row. Sample IDs must match PLINK .fam file exactly.
  Each locus requires two columns: <locus>_h1 and <locus>_h2

  Example:
    Sample,KIR3DL1_h1,KIR3DL1_h2
    HG00096,KIR3DL1*001,KIR3DL1*002
    HG00097,KIR3DL1*005,KIR3DL1*015
    HG00099,KIR3DL1*0000,KIR3DL1*020   ← null allele = *0000

  ${YELLOW}Note: rows with *new or *unresolved alleles are automatically excluded${RESET}

${BOLD}${WHITE}EXAMPLES${RESET}
  # Basic training
  pong2 train -i data/ref_chr19 -k data/kir_calls.csv -o models/ -l KIR3DL1 -a hg19

  # With custom parameters
  pong2 train -i ref_chr19 -k kir_calls.csv -o models/ -l KIR3DL1 -a hg38 \\
              -t 20 --nclassifier 200 --split 0.8 --kirmaf 0.005

${BOLD}${WHITE}OUTPUT FILES${RESET}
  <locus>_model.RData      Trained prediction model (main output)
  <locus>_test.RData       Test genotypes (only when --split < 1)
  <locus>_split.RData      Train/test split object (only when --split < 1)

${BOLD}${WHITE}TROUBLESHOOTING${RESET}
  • No matching samples       → Sample IDs in --kfile must exactly match .fam file
  • Insufficient samples      → Need ≥ 10 overlapping samples between KIR and PLINK files
  • No SNPs in region         → Check --assembly and --region coordinates
  • Slow training             → Increase --threads or reduce --nclassifier
  • Low model accuracy        → Increase sample size or adjust --kirmaf threshold
HELP
}

# ────────────────────────────────────────────────
# Dispatcher
# ────────────────────────────────────────────────
case "${1:-}" in
    "impute")
        show_impute_help
        ;;
    "train")
        show_train_help
        ;;
    "" | "-h" | "--help" | "help")
        show_main_help
        ;;
    "version" | "--version" | "-v")
        echo "PONG2 version 1.0.0"
        ;;
    *)
        echo -e "${RED}Error: Unknown command or option '$1'${RESET}" >&2
        echo -e "Use: pong2 --help" >&2
        exit 1
        ;;
esac

exit 0
