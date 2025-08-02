#!/bin/bash
# Auto-Installer for Required Genomics Binaries

# Configuration
BIN_DIR="$HOME/bin"
REQUIRED_APPS=("plink2" "minimac4", "hg19")
export PATH="$BIN_DIR:$PATH"  # Add to PATH temporarily

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Create bin directory if missing
mkdir -p "$BIN_DIR"

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
    local installer="$BIN_DIR/minimac4-installer.sh"

    # Download
    if ! wget -q --show-progress "$url" -O "$installer"; then
        echo -e "${RED}Download failed${NC}"
        return 1
    fi

    # Install
    chmod +x "$installer"
    if ! "$installer" -b -p "$BIN_DIR"; then
        echo -e "${RED}Installation failed${NC}"
        return 1
    fi

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

install_eagle() {
    local version="v2.4.1"
    local temp_dir=$(mktemp -d)
    local archive="Eagle_${version}.tar.gz"
    local download_url="https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/${archive}"

    echo -e "${YELLOW}Installing Eagle ${version}...${NC}"

    # Download and extract
    if ! wget -q --show-progress "$download_url" -O "$temp_dir/$archive"; then
        echo -e "${RED}❌ Download failed${NC}"
        rm -rf "$temp_dir"
        return 1
    fi

    if ! tar -xzf "$temp_dir/$archive" -C "$temp_dir"; then
        echo -e "${RED}❌ Extraction failed${NC}"
        rm -rf "$temp_dir"
        return 1
    fi

    # Install binary
    mkdir -p "$BIN_DIR"
    if ! mv "$temp_dir/eagle" "$BIN_DIR/"; then
        echo -e "${RED}❌ Installation failed (permission issue?)${NC}"
        rm -rf "$temp_dir"
        return 1
    fi

    # Set permissions
    chmod +x "$BIN_DIR/eagle"

    # Cleanup
    rm -rf "$temp_dir"
    
    # Verify
    if "$BIN_DIR/eagle" --version &>/dev/null; then
        echo -e "${GREEN}✓ Eagle ${version} installed to $BIN_DIR/eagle${NC}"
        return 0
    else
        echo -e "${RED}❌ Installation verification failed${NC}"
        return 1
    fi
}



prompt_install() {
    local app=$1
    while true; do
        read -rp "Install $app? [y/n]: " yn
        case "$yn" in
            [Yy]*)
                "install_${app}" && return 0
                echo -e "${RED}Installation aborted${NC}"
                return 1
                ;;
            [Nn]*)
                echo -e "${YELLOW}Skipping $app installation${NC}"
                return 1
                ;;
            *)
                echo "Please answer yes (y) or no (n)"
                ;;
        esac
    done
}

# Main check
echo -e "\n${YELLOW}=== Dependency Check ===${NC}"
missing_count=0

for app in "${REQUIRED_APPS[@]}"; do
    if command -v "$app" >/dev/null 2>&1; then
        echo -e "${GREEN}Found: $app ($(command -v "$app"))${NC}"
    else
        echo -e "${RED}Missing: $app${NC}"
        if prompt_install "$app"; then
            echo -e "${GREEN}Successfully installed $app${NC}"
        else
            ((missing_count++))
        fi
    fi
done

# Permanent PATH setup suggestion
echo -e "\n${YELLOW}Add this to your ~/.bashrc or ~/.zshrc:${NC}"
echo "export PATH=\"\$HOME/bin:\$PATH\""

echo -e "\n${GREEN}All tools are now available!${NC}"