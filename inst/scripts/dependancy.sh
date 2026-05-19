#!/usr/bin/env bash
# =============================================================================
# PONG2 Dependency Checker & Auto-Installer
# Checks and installs: plink2, minimac4, and 1KGP reference panel
# =============================================================================
set -euo pipefail

# =============================================================================
# ARGUMENTS
# =============================================================================
OUTPUT_DIR="$1"
ASSEMBLY="$2"
TOOL="$3"

# =============================================================================
# CONFIGURATION
# =============================================================================
BIN_DIR="$HOME/bin"
mkdir -p "$BIN_DIR"
mkdir -p "$OUTPUT_DIR/tmp"
export PATH="$BIN_DIR:$PATH"

PONG2_root=$(Rscript -e 'cat(system.file(package="PONG2"))' 2>/dev/null)

# ── Reference file names per assembly ────────────────────────────────────────
HG19_REF="ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes"
HG38_REF="1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel"

if [[ "$ASSEMBLY" == "hg19" ]]; then
  ref_file="$HG19_REF"
  ref_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${HG19_REF}.vcf.gz"
  ref_tbi_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${HG19_REF}.vcf.gz.tbi"
else
  ref_file="$HG38_REF"
  ref_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/${HG38_REF}.vcf.gz"
  ref_tbi_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/${HG38_REF}.vcf.gz.tbi"
fi

ref_dir="$PONG2_root/extdata/$ASSEMBLY"
ref_msav="$ref_dir/$ref_file.msav"

# ── Required tools ────────────────────────────────────────────────────────────
if [[ "$TOOL" == "minimac4" ]]; then
  REQUIRED_APPS=("minimac4" "hg")
else
  REQUIRED_APPS=("$TOOL")
fi

# =============================================================================
# COLORS
# readable on both white and dark terminal backgrounds
# =============================================================================
RED='\033[0;31m'      # errors
GREEN='\033[0;32m'    # success
CYAN='\033[0;36m'     # installing / downloading (replaces faded yellow)
BLUE='\033[0;34m'     # info / skipping
BOLD='\033[1m'        # section headers / emphasis
NC='\033[0m'

# =============================================================================
# INSTALL FUNCTIONS
# =============================================================================
install_plink2() {
  echo -e "${CYAN}Installing PLINK2...${NC}"

  case "$(uname -s)-$(uname -m)" in
    Linux-x86_64)
      if grep -q avx2 /proc/cpuinfo 2>/dev/null; then
        url="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20250707.zip"
      else
        url="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250707.zip"
      fi
      ;;
    Darwin-arm64)  url="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_arm64_20250707.zip" ;;
    Darwin-x86_64) url="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_x86_64_20250707.zip" ;;
    CYGWIN*|MINGW*|MSYS*)
      url="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_win64_20250707.zip" ;;
    *)
      echo -e "${RED}❌  Unsupported platform: $(uname -s)-$(uname -m)${NC}"
      echo -e "${BLUE}Manual download: https://www.cog-genomics.org/plink/2.0/${NC}"
      return 1
      ;;
  esac

  echo -e "${BLUE}Downloading: $url${NC}"
  if ! wget -q --show-progress "$url" -O "$OUTPUT_DIR/tmp/plink2.zip"; then
    echo -e "${RED}❌  Download failed — check network or download manually:${NC}"
    echo -e "${BLUE}https://www.cog-genomics.org/plink/2.0/${NC}"
    return 1
  fi

  if ! unzip -qo "$OUTPUT_DIR/tmp/plink2.zip" -d "$BIN_DIR"; then
    echo -e "${RED}❌  Extraction failed (corrupted download)${NC}"
    return 1
  fi

  if [[ ! "$(uname -s)" =~ CYGWIN|MINGW|MSYS ]]; then
    chmod +x "$BIN_DIR"/plink2*
  fi

  rm -f "$OUTPUT_DIR/tmp/plink2.zip"
  echo -e "${GREEN}✅  PLINK2 installed to $BIN_DIR/plink2${NC}"
}

install_minimac4() {
  if [[ "$(uname -s)" != "Linux" ]] || [[ "$(uname -m)" != "x86_64" ]]; then
    echo -e "${RED}❌  Minimac4 binary only supports Linux x86_64${NC}"
    echo -e "${BLUE}Build from source: https://github.com/statgen/Minimac4/releases/tag/v4.1.6${NC}"
    return 1
  fi

  echo -e "${CYAN}Installing Minimac4 v4.1.6...${NC}"
  local url="https://github.com/statgen/Minimac4/releases/download/v4.1.6/minimac4-4.1.6-Linux-x86_64.sh"
  local installer="$OUTPUT_DIR/tmp/minimac4-installer.sh"

  if ! wget -q --show-progress "$url" -O "$installer"; then
    echo -e "${RED}❌  Download failed${NC}"
    return 1
  fi

  chmod +x "$installer"
  if ! "$installer" --skip-license --prefix="$OUTPUT_DIR/tmp"; then
    echo -e "${RED}❌  Installation failed${NC}"
    return 1
  fi

  mv "$OUTPUT_DIR/tmp/bin/minimac4" "$BIN_DIR/"
  rm -f "$installer"

  if "$BIN_DIR/minimac4" --version &>/dev/null; then
    echo -e "${GREEN}✅  Minimac4 installed to $BIN_DIR/minimac4${NC}"
    return 0
  else
    echo -e "${RED}❌  Minimac4 installation verification failed${NC}"
    return 1
  fi
}

install_hg() {
  echo -e "${CYAN}Downloading ${ASSEMBLY} reference panel...${NC}"

  local vcf_out="$OUTPUT_DIR/tmp/${ref_file}.vcf.gz"
  local tbi_out="$OUTPUT_DIR/tmp/${ref_file}.vcf.gz.tbi"

  echo -e "${BLUE}Downloading VCF: $ref_url${NC}"
  if ! wget -q --show-progress "$ref_url" -O "$vcf_out"; then
    echo -e "${RED}❌  VCF download failed${NC}"
    return 1
  fi

  echo -e "${BLUE}Downloading index: $ref_tbi_url${NC}"
  if ! wget -q --show-progress "$ref_tbi_url" -O "$tbi_out"; then
    echo -e "${RED}❌  Index download failed${NC}"
    return 1
  fi

  echo -e "${CYAN}Compressing reference to msav format (this may take several minutes)...${NC}"
  mkdir -p "$ref_dir"

  if ! "$BIN_DIR/minimac4" --compress-reference "$vcf_out" > "$OUTPUT_DIR/tmp/${ref_file}.msav"; then
    echo -e "${RED}❌  minimac4 --compress-reference failed${NC}"
    return 1
  fi

  if ! mv "$OUTPUT_DIR/tmp/${ref_file}.msav" "$ref_msav"; then
    echo -e "${RED}❌  Failed to move msav to $ref_dir${NC}"
    return 1
  fi

  rm -f "$vcf_out" "$tbi_out"
  echo -e "${GREEN}✅  Reference panel saved to $ref_msav${NC}"
}

# =============================================================================
# PROMPT INSTALL
# =============================================================================
prompt_install() {
  local app="$1"
  local label

  if [[ "$app" == "hg" ]]; then
    label="Download $ASSEMBLY reference panel"
  else
    label="Install $app"
  fi

  while true; do
    read -rp "${label}? [y/n]: " yn
    case "$yn" in
      [Yy]*)
        if [[ "$app" == "hg" ]]; then
          install_hg
        else
          install_"${app}"
        fi
        return $?
        ;;
      [Nn]*)
        echo -e "${BLUE}Skipping $label${NC}"
        return 1
        ;;
      *)
        echo "Please answer yes (y) or no (n)"
        ;;
    esac
  done
}

# =============================================================================
# MAIN DEPENDENCY CHECK
# =============================================================================
echo -e "\n${BOLD}=== Dependency Check ===${NC}"
missing_count=0

for app in "${REQUIRED_APPS[@]}"; do

  # ── Check hg reference panel ───────────────────────────────────────────────
  if [[ "$app" == "hg" ]]; then
    if [[ -f "$ref_msav" ]]; then
      echo -e "${GREEN}✅  Found: $ASSEMBLY reference ($ref_msav)${NC}"
    else
      echo -e "${RED}Missing: $ASSEMBLY reference panel${NC}"
      if ! prompt_install "hg"; then
        missing_count=$((missing_count + 1))
        echo -e "${RED}❌  $ASSEMBLY reference panel required for --fill-missing${NC}"
        exit 1
      fi
    fi

  # ── Check binary tools ─────────────────────────────────────────────────────
  else
    if command -v "$app" >/dev/null 2>&1; then
      echo -e "${GREEN}✅  Found: $app ($(command -v "$app"))${NC}"
    else
      echo -e "${RED}Missing: $app${NC}"
      if ! prompt_install "$app"; then
        missing_count=$((missing_count + 1))
        echo -e "${RED}❌  $app is required but not installed${NC}"
        exit 1
      fi
    fi
  fi

done

# ── Final PATH export ─────────────────────────────────────────────────────────
export PATH="$HOME/bin:$PATH"

if [[ $missing_count -eq 0 ]]; then
  echo -e "\n${GREEN}✅  All dependencies satisfied${NC}"
else
  echo -e "\n${RED}❌  $missing_count dependency/dependencies missing${NC}"
  exit 1
fi
