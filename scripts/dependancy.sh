#!/bin/bash
# Auto-Installer for Required Genomics Binaries

OUTPUT_DIR="$1"
ASSEMBLY="$2"
TOOL=$3
# Configuration
BIN_DIR="$HOME/bin"

if [ "$TOOL" == "minimac4" ]; then
  REQUIRED_APPS=("$TOOL" "hg")
else
  REQUIRED_APPS=("$TOOL")
fi

export PATH="$BIN_DIR:$PATH"  # Add to PATH temporarily

hg=("ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes", "1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel")

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Create bin directory if missing
mkdir -p "$BIN_DIR"
mkdir -p "$OUTPUT_DIR/tmp"
PONG2_root=$(Rscript -e 'cat(system.file(package="PONG2"))' 2>/dev/null)

if [ $ASSEMBLY = "hg19" ]; then
  ref_file="ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes"
else
  ref_file="1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel"
fi


# Version requirements configuration
declare -A VERSION_REQUIREMENTS=(
    ["minimac4"]="4.1.6"
    ["plink2"]="2.0"
    ["bcftools"]="1.10"
)

# Function to get current version
get_version() {
    local app="$1"
    case "$app" in
        "minimac4")
            minimac4 --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1
            ;;
        "plink2")
            plink2 --version 2>/dev/null | grep -oE 'v[0-9]+\.[0-9]+' | sed 's/v//' | head -1
            ;;
        "bcftools")
            bcftools --version 2>/dev/null | head -1 | grep -oE '[0-9]+\.[0-9]+' | head -1
            ;;
        *)
            echo "unknown"
            ;;
    esac
}

# Function to compare versions numerically (more robust)
version_compare() {
    local current="$1"
    local required="$2"

    # Using sort -V for proper version comparison
    if printf "%s\n%s" "$current" "$required" | sort -V -C 2>/dev/null; then
        return 0  # current >= required
    else
        return 1  # current < required
    fi
}

# Function to check version match
check_version_match() {
    local app="$1"
    local required_version="${VERSION_REQUIREMENTS[$app]}"

    # If no version requirement, return true
    [[ -z "$required_version" ]] && return 0

    local current_version=$(get_version "$app")

    # If version check failed, return false
    [[ "$current_version" == "unknown" ]] && return 1

    # Compare versions (simple string comparison - adjust if needed)
    [[ "$current_version" == "$required_version" ]]
}



install_plink2() {
    echo -e "${YELLOW}Installing PLINK2...${NC}"

    # Current working URLs (verified July 2024)
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64)
            if grep -q avx2 /proc/cpuinfo; then
                url="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20250707.zip"
            else
                url="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250707.zip"
            fi
            ;;
        Darwin-arm64) url="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_arm64_20250707.zip" ;;
        Darwin-x86_64) url="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_x86_64_20250707.zip" ;;
        CYGWIN*-*|MINGW*-*|MSYS*-*)  # Windows (Cygwin/Git Bash)
            url="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_win64_20250707.zip"
            ;;
        *)
            echo -e "${RED}❌ Unsupported platform: $(uname -s)-$(uname -m)${NC}"
            return 1
            ;;
    esac

    # Download with wget
    echo -e "Downloading: ${BLUE}$url${NC}"
    if ! wget -q --show-progress "$url" -O plink2.zip; then
        echo -e "${RED}❌ Download failed (404 or network issue)${NC}"
        echo -e "Try manual download from:"
        echo -e "https://www.cog-genomics.org/plink/2.0/"
        return 1
    fi

    # Extract
    if ! unzip -qo plink2.zip -d "$BIN_DIR"; then
        echo -e "${RED}❌ Extraction failed (corrupted download)${NC}"
        return 1
    fi

    # Set permissions (except Windows)
    if [[ ! "$(uname -s)" =~ CYGWIN|MINGW|MSYS ]]; then
        chmod +x "$BIN_DIR"/plink2*
    fi

    rm plink2.zip
    echo -e "${GREEN}✓ PLINK2 installed to $BIN_DIR/${NC}"
    echo -e "Add to PATH: ${YELLOW}export PATH=\"\$PATH:$BIN_DIR\"${NC}"
}

install_minimac4() {
    # Only support Linux x86_64
    if [[ "$(uname -s)" != "Linux" ]] || [[ "$(uname -m)" != "x86_64" ]]; then
        echo -e "${RED}Error: Minimac4 only provides Linux x86_64 binaries${NC}"
        echo -e "Try building from the source:"
        echo -e "https://github.com/statgen/Minimac4/releases/tag/v4.1.6"
        return 1
    fi

    echo -e "${YELLOW}Installing Minimac4 for Linux...${NC}"
    local url="https://github.com/statgen/Minimac4/releases/download/v4.1.6/minimac4-4.1.6-Linux-x86_64.sh"
    local installer="$OUTPUT_DIR/tmp/minimac4-installer.sh"

    # Download
    if ! wget -q --show-progress "$url" -O "$installer"; then
        echo -e "${RED}Download failed${NC}"
        return 1
    fi

    # Install
    chmod +x "$installer"
    if ! "$installer" --skip-license --prefix="$OUTPUT_DIR/tmp"; then
        echo -e "${RED}Installation failed${NC}"
        return 1
    fi

    #move minimac to  $BIN_DIR
    mv "$OUTPUT_DIR/tmp/bin/minimac4" $BIN_DIR/

    # Verify
    if "$BIN_DIR/minimac4" --version &>/dev/null; then
        echo -e "${GREEN}✓ Minimac4 installed to $BIN_DIR/minimac4${NC}"
        rm "$installer"
        return 0
    else
        echo -e "${RED}Installation verification failed${NC}"
        return 1
    fi
}

install_ref(){
    # Check if ASSEMBLY variable equals "hg19"
    if [ "$ASSEMBLY" = "hg19" ]; then
        local ref="$PONG2_root/extdata/hg19/"
        local installer1="$OUTPUT_DIR/tmp/${hg[0]}.vcf.gz"
        local installer2="$OUTPUT_DIR/tmp/${hg[0]}.vcf.gz.tbi"

        echo -e "${YELLOW}Downloading hg19 reference...${NC}"
        local url1="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${hg[0]}.vcf.gz"
        local url2="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${hg[0]}.vcf.gz.tbi"

    else
        local ref="$PONG2_root/extdata/hg38/"
        local installer1="$OUTPUT_DIR/tmp/${hg[1]}.vcf.gz"
        local installer2="$OUTPUT_DIR/tmp/${hg[1]}.vcf.gz.tbi"

        echo -e "${YELLOW}Downloading hg38 reference...${NC}"
        local url1="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/${hg[1]}.vcf.gz"
        local url2="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/${hg[1]}.vcf.gz.tbi"

    fi


    # Download
    if ! wget -q --show-progress "$url1" -O "$installer1"; then
        echo -e "${RED}Download failed${NC}"
        return 1
    fi

    if ! wget -q --show-progress "$url2" -O "$installer2"; then
        echo -e "${RED}Download failed${NC}"
        return 1
    fi

    #minimac4 --compress-reference reference.{sav,bcf,vcf.gz} > reference.msav
    echo -e "${BLUE}Compressing reference data (this may take several minutes)...${NC}"
    if ! "$BIN_DIR/minimac4" "--compress-reference" "$installer1" > "$OUTPUT_DIR/tmp/$ref_file.msav"; then
      echo -e "${RED}minimac4 failed for --compress-reference ${NC}"
      return 1
    fi

    if ! mv "$OUTPUT_DIR/tmp/$ref_file.msav" "$ref"; then
      return 1
    fi
}

prompt_install() {
    local app=$1
    while true; do
        if [[ "$app" == "hg" ]]; then
            Install="Download $ASSEMBLY"
            Installation="Downloading"
        else
            Install="Install $app"
            Installation="Installation"
        fi
        read -rp "$Install? [y/n]: " yn
        case "$yn" in
            [Yy]*)
                if [[ "$app" == "hg" ]]; then
                  "install_ref"
                else
                    "install_${app}"
                fi

                # Check if installation succeeded
                if [[ $? -eq 0 ]]; then
                    return 0
                else
                    echo -e "${RED}$Installation failed for $app${NC}"
                    return 1
                fi
                ;;
            [Nn]*)
                echo -e "${YELLOW}Skipping $app $Installation${NC}"
                return 1
                ;;
            *)
                echo "Please answer yes (y) or no (n)"
                ;;
        esac
    done
}

# Main check
echo -e "\n${GREEN}=== Dependency Check ===${NC}"
missing_count=0

for app in "${REQUIRED_APPS[@]}"; do
    if [[ "$app" == "plink2" || "$app" == "minimac4" ]] && command -v "$app" >/dev/null 2>&1 && check_version_match "$app"; then
        version=$(get_version "$app")
        required="${VERSION_REQUIREMENTS[$app]:-none}"
        echo -e "${GREEN}✓ Found: $app, meets requirement: v$required ($(command -v "$app"))${NC}"

    elif [[ "$app" == "hg" ]] && [[ -f "$PONG2_root/extdata/$ASSEMBLY/$ref_file.msav" ]]; then
      echo -e "${GREEN}✓ Found: compressed genome build ($PONG2_root/extdata/$ASSEMBLY/$ref_file.msav)${NC}"

    # App not found?
    else
       if [[ "$app" != "${hg[0]}" && "$app" != "${hg[1]}" ]]; then
          appv="${VERSION_REQUIREMENTS[$app]}"
          echo -e "${RED}Missing: $app  v$appv ${NC}"
          Install="install $app v$appv"
       else
         echo -e "${GREEN} Missing genome build specification ($ASSEMBLY)${NC}"
         Install="download $ASSEMBLY"
       fi
       if prompt_install "$app"; then
            echo -e "${GREEN}Successfully $Install${NC}"
        else
            ((missing_count++))
            echo $missing_count
            exit 1
        fi
    fi
done

# Permanent PATH setup suggestion
#echo -e "\n${YELLOW}Add this to your ~/.bashrc or ~/.zshrc:${NC}"
export PATH=$HOME/bin:$PATH

#echo -e "\n${GREEN}All tools are now available!${NC}"
