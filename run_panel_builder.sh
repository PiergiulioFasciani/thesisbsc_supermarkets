#!/bin/bash

# Colors
BLUE='\033[0;34m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m' # No Color

echo ""
echo "================================================================================"
echo -e "${BOLD}${BLUE}OSM PANEL BUILDER${NC}"
echo "================================================================================"
echo ""

# Check for config file
CONFIG_FILE="config.yml"
if [ ! -f "$CONFIG_FILE" ]; then
    echo -e "${RED}✗ Config file not found: $CONFIG_FILE${NC}"
    exit 1
fi

echo -e "${YELLOW}→${NC} Using config: ${BOLD}$CONFIG_FILE${NC}"

# Get the sample directory (most recent or specified)
if [ -z "$1" ]; then
    echo ""
    echo -e "${YELLOW}→${NC} Config will auto-detect most recent sample"
else
    SAMPLE_DIR="$1"
    if [ ! -d "$SAMPLE_DIR" ]; then
        echo -e "${RED}✗ Directory not found: $SAMPLE_DIR${NC}"
        exit 1
    fi
    INPUT_FILE="${SAMPLE_DIR}/establishments.csv"
    if [ ! -f "$INPUT_FILE" ]; then
        echo -e "${RED}✗ Establishments file not found: $INPUT_FILE${NC}"
        exit 1
    fi
    echo ""
    echo -e "${YELLOW}→${NC} Using specified sample: ${BOLD}$INPUT_FILE${NC}"
fi

echo ""
echo "--------------------------------------------------------------------------------"
echo -e "${BOLD}STEP 1: Activate Environment${NC}"
echo "--------------------------------------------------------------------------------"

# Activate conda environment
if [ -d "$HOME/miniconda3" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate osm_sampler_env
    echo -e "${GREEN}✓${NC} Environment activated: osm_sampler_env"
else
    echo -e "${RED}✗ Miniconda not found${NC}"
    exit 1
fi

echo ""
echo "--------------------------------------------------------------------------------"
echo -e "${BOLD}STEP 2: Generate Panel${NC}"
echo "--------------------------------------------------------------------------------"
echo ""

# Run panel builder
if [ -n "$INPUT_FILE" ]; then
    python scripts/panel_builder.py --config "$CONFIG_FILE" --input "$INPUT_FILE"
else
    python scripts/panel_builder.py --config "$CONFIG_FILE"
fi

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "================================================================================"
    echo -e "${BOLD}${GREEN}✓ SUCCESS${NC}"
    echo "================================================================================"
    echo ""
    echo -e "${GREEN}→${NC} Panel dataset generated successfully!"
    echo -e "${GREEN}→${NC} Check data/output/panels/ for the results"
    echo ""
else
    echo ""
    echo "================================================================================"
    echo -e "${BOLD}${RED}✗ FAILED${NC}"
    echo "================================================================================"
    echo ""
    echo -e "${RED}→${NC} Panel generation failed with exit code $EXIT_CODE"
    echo -e "${RED}→${NC} Check the error messages above"
    echo ""
fi

exit $EXIT_CODE
