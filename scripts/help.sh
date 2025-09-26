#!/bin/bash
# PONG2 Help System
# KIR Imputation and Analysis Toolkit

# Detect if running as R package or CLI
if [[ "$0" == *"inst/scripts/help.sh" ]] || [[ -n "$R_PACKAGE_DIR" ]]; then
    MODE="R"
else
    MODE="CLI"
fi

show_help() {
    # Color definitions
    local BOLD=$(tput bold)
    local GREEN=$(tput setaf 2)
    local BLUE=$(tput setaf 4)
    local CYAN=$(tput setaf 6)
    local WHITE=$(tput setaf 7)
    local YELLOW=$(tput setaf 3)
    local RED=$(tput setaf 1)
    local RESET=$(tput sgr0)

    cat << EOF
${BOLD}${GREEN}PONG2 - KIR Genotype Imputation and Analysis Toolkit${RESET}
${BLUE}Version: ${WHITE}2.1.0${RESET}
${BLUE}Description: ${RESET}Advanced KIR imputation pipeline with quality control and visualization

${BOLD}${GREEN}USAGE:${RESET}
  pong2 impute [OPTIONS] -i <input> -o <output> -l <locus> -a <assembly>

${BOLD}${GREEN}REQUIRED OPTIONS:${RESET}
  ${CYAN}-i, --input FILE${RESET}       Input file plink/vcf path (phased if using minimac4)
  ${CYAN}-o, --output DIR${RESET}       Output directory path
  ${CYAN}-l, --locus GENE${RESET}       Target KIR locus (e.g., KIR2DL1, KIR3DL1, KIR2DS4)
  ${CYAN}-a, --assembly ASSEMBLY${RESET} Reference genome assembly (hg19, hg38)

${BOLD}${GREEN}OTHER OPTIONS:${RESET}
  ${CYAN}-t, --threads NUM${RESET}      Number of CPU threads [default: 4]
  ${CYAN}-f, --force${RESET}            force the imputation
  ${CYAN}--filter VALUE${RESET}         Filter threshold (0.01, 0.005) [default: 0.01]
  ${CYAN}--fill-missing${RESET}         Fill missing genotypes in the input file with minimac4
  ${CYAN}-h, --help${RESET}             Show this help message

${BOLD}${GREEN}EXAMPLES:${RESET}
  ${WHITE}Basic imputation:${RESET}
    pong2 impute -i chr19 \\
              -o ~/out -l KIR2DL1 -a hg19

  ${WHITE}With advanced options:${RESET}
    pong2 impute -i chr19_phased -o /out \\
              -l KIR2DL1 -a hg38 --threads 8 --fill-missing --filter 0.005

  ${WHITE}Force overwrite:${RESET}
    pong2 impute --input ./imputed_data --output ./new_results \\
              --locus KIR3DL1 --assembly hg19 --force

${BOLD}${GREEN}SUPPORTED LOCI:${RESET}
 KIR2DL1, KIR2DL3, KIR2DL4, KIR3DL1, KIR3DL2, KIR2DS1, KIR2DS2, KIR2DS4, KIR2DS5

${BOLD}${YELLOW}INPUT REQUIREMENTS:${RESET}
  • Phased imputation data in PLINK format
  • Chromosome 19 KIR region data
  • Properly formatted sample IDs and coordinates

${BOLD}${RED}TROUBLESHOOTING:${RESET}
  • Ensure input files are properly phased and formatted
  • Verify reference assembly matches your data
  • Check that output directory is writable
  • Use --force to overwrite existing results

${BOLD}${BLUE}DOCUMENTATION:${RESET}
  Full documentation: https://github.com/NormanLabUCD/PONG2
  Issue tracking: https://github.com/NormanLabUCD/PONG2/issues

${BOLD}${BLUE}CITATION:${RESET}
  Please cite our publication when using PONG2 in your research.

  Sadeeq, S. A., Leaton, L., Castelli, E., & Norman, P. (2025).PONG 2.0:
  Allele imputation for the killer cell immunoglobulin-like receptors. Human Immunology, 86, 111488.
EOF
}

show_impute_help() {
    local BOLD=$(tput bold)
    local GREEN=$(tput setaf 2)
    local CYAN=$(tput setaf 6)
    local WHITE=$(tput setaf 7)
    local RESET=$(tput sgr0)

    cat << EOF

${BOLD}${GREEN}Imputation Command Details${RESET}

${BOLD}${WHITE}Workflow:${RESET}
  1. Load phased imputation data for specified locus
  2. Apply quality filters (--filter parameter)
  3. Perform KIR genotype imputation
  4. Generate output files with optional missing data filling

${BOLD}${WHITE}Input Format:${RESET}
  Expected: Phased imputation output for chromosome 19 KIR region
  Required files:
    • imputed_chr19_kir_phased (or similar)
    • Proper sample alignment with reference

${BOLD}${WHITE}Output Files:${RESET}
  • <output>/<locus>_<assembly>.txt - Main imputation results
  • <output>/quality_metrics.csv - Quality control metrics
  • <output>/imputation_stats.json - Statistical summary

${BOLD}${WHITE}Quality Filtering:${RESET}
  The --filter parameter controls imputation quality threshold:
  • 0.01: High quality (more stringent)
  • 0.005: Permissive (for low-quality data)

${BOLD}${WHITE}Missing Data Handling:${RESET}
  Use --fill-missing to impute missing genotypes rather than excluding them
  This increases completeness but may reduce accuracy for low-quality samples

${BOLD}${WHITE}Performance Tips:${RESET}
  • Use --threads to parallelize for large datasets
  • Ensure sufficient memory for large sample sizes
  • --force avoids prompt for overwriting existing files
EOF
}

# Main help dispatcher
case "${1:-}" in
    "impute")
        show_help
        show_impute_help
        ;;
    "qc"|"quality"|"qualitycheck")
        echo -e "\nQuality control help not yet implemented"
        ;;
    "plot"|"visualize")
        echo -e "\nVisualization help not yet implemented"
        ;;
    "version")
        echo "PONG2 v2.1.0"
        ;;
    "-h"|"--help"|"help"|"")
        show_help
        ;;
    *)
        echo -e "${RED}Error: Unknown command or option '$1'${NC}" >&2
        echo -e "Use 'pong2 --help' for usage information" >&2
        exit 1
        ;;
esac

exit 0
