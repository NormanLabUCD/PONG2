#!/bin/bash
# PONG2 Imputation Module (PLINK version)
# Usage: ./impute.sh --input <prefix> --output <dir> [OPTIONS]

set -euo pipefail  # Exit on error/undefined vars

# --------------------------
# Configuration
# --------------------------


# --------------------------
# Color Definitions
# --------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

#  Smart Path Detection ---------------------------------------------------
get_package_root() {
    # Try package installation path first
    local pkg_path=$(Rscript -e 'cat(system.file(package="PONG"))' 2>/dev/null)
    
    if [[ -z "$pkg_path" ]]; then
        # Development mode - get path relative to this script
        local script_dir=$(dirname "$(realpath "$0")")
        echo "$(realpath "$script_dir/../..")"  # Move up from inst/scripts/
    else
        echo "$pkg_path"
    fi
}

PACKAGE_ROOT=$(get_package_root)


# --------------------------
# Functions
# --------------------------
die() {
  echo -e "${RED}ERROR: $*${NC}" >&2
  exit 1
}

validatePlink() {
    local prefix="$1"
    local output="$2"

    echo -e "${GREEN}Validating PLINK files...${NC}"
     [[ -f "${prefix}.bed" && -f "${prefix}.fam" && -f "${prefix}.bim" ]] || 
    die "Missing PLINK files for prefix: $prefix"

    # Create output directory if it doesn't exist
    mkdir -p "$output" || die "Cannot create output directory"
    mkdir -p "$output/temp" || die "Cannot create output directory"

    awk '$1 == 19 && $4 >= 55235681 && $4 <= 55378697 {print}' "$prefix.bim" > "$output/temp/KIR_region.txt" || die "No SNPs found in the KIR chromosome 19 region"

    awk 'NR==FNR {ref[$1":"$2]=1; next} 
     {pos_key=$1":"$4} 
     pos_key in ref' "$PACKAGE_ROOT/data/hg19/kirhg19_pos.txt" "$output/temp/KIR_region.txt" > "$output/temp/matched_snps.txt"
    
    matched_count=$(wc -l < "$output/temp/matched_snps.txt")
    ref_count=$(wc -l < "$PACKAGE_ROOT/data/hg19/kirhg19_pos.txt")
    match_percent=$(awk -v matched="$matched_count" -v total="$ref_count" 'BEGIN {printf "%.1f", (matched/total)*100}')
}
