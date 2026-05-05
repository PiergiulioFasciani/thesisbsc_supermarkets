#!/bin/bash
# Extract Lombardy region from Italy OSM history file
# One-time extraction to create a smaller, focused dataset for Milan study

set -e

# Colors
GREEN='\033[0;32m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m'

echo ""
echo "================================================================================"
echo -e "${BOLD}${CYAN}OSM HISTORY EXTRACTION - LOMBARDY REGION${NC}"
echo "================================================================================"
echo ""

# Configuration
INPUT_FILE="data/input/italy_06022025.osh.pbf"
OUTPUT_FILE="data/input/lombardy.osh.pbf"

# Lombardy bounding box (approximate)
# Covers Milan metropolitan area and surrounding Lombardy region
MIN_LON="8.4"
MIN_LAT="44.8"
MAX_LON="10.6"
MAX_LAT="46.1"

echo -e "${CYAN}→${NC} Input file: ${BOLD}${INPUT_FILE}${NC}"
echo -e "${CYAN}→${NC} Output file: ${BOLD}${OUTPUT_FILE}${NC}"
echo -e "${CYAN}→${NC} Bounding box: [${MIN_LON}, ${MIN_LAT}, ${MAX_LON}, ${MAX_LAT}]"
echo ""

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${YELLOW}WARNING:${NC} Input file not found: ${INPUT_FILE}"
    echo "Please ensure the Italy OSM history file is in data/input/"
    exit 1
fi

# Get input file size
INPUT_SIZE=$(du -h "$INPUT_FILE" | cut -f1)
echo -e "${CYAN}→${NC} Input file size: ${BOLD}${INPUT_SIZE}${NC}"
echo ""

# Check if osmium-tool is installed
echo -ne "Checking for osmium-tool... "
if ! command -v osmium &> /dev/null; then
    echo -e "${YELLOW}NOT FOUND${NC}"
    echo ""
    echo -e "${YELLOW}osmium-tool is not installed.${NC}"
    echo ""
    echo "Install it with:"
    echo -e "  ${CYAN}brew install osmium-tool${NC}  (macOS)"
    echo -e "  ${CYAN}sudo apt install osmium-tool${NC}  (Ubuntu/Debian)"
    echo -e "  ${CYAN}conda install -c conda-forge osmium-tool${NC}  (Conda)"
    echo ""
    exit 1
fi
echo -e "${GREEN}✓${NC}"
echo ""

echo "================================================================================"
echo -e "${BOLD}EXTRACTING LOMBARDY REGION${NC}"
echo "================================================================================"
echo ""
echo "This may take 30-60 minutes for a 3.8GB history file..."
echo ""

# For OSM history files, we need to use a different approach
# osmium extract doesn't work with history files directly
# We'll use osmium tags-filter combined with osmium getid

echo "Step 1: Filtering by geographic area (this will take a while)..."
echo ""

# Extract using osmium extract with history option
osmium extract \
    --bbox ${MIN_LON},${MIN_LAT},${MAX_LON},${MAX_LAT} \
    --with-history \
    --strategy complete_ways \
    --set-bounds \
    "${INPUT_FILE}" \
    -o "${OUTPUT_FILE}" \
    --overwrite

echo ""
echo "================================================================================"
echo -e "${GREEN}${BOLD}✓ EXTRACTION COMPLETE${NC}"
echo "================================================================================"
echo ""

# Get output file size
OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
echo -e "${CYAN}→${NC} Output file: ${BOLD}${OUTPUT_FILE}${NC}"
echo -e "${CYAN}→${NC} Output file size: ${BOLD}${OUTPUT_SIZE}${NC}"
echo ""

# Calculate size reduction
INPUT_BYTES=$(stat -f%z "$INPUT_FILE" 2>/dev/null || stat -c%s "$INPUT_FILE" 2>/dev/null)
OUTPUT_BYTES=$(stat -f%z "$OUTPUT_FILE" 2>/dev/null || stat -c%s "$OUTPUT_FILE" 2>/dev/null)
REDUCTION=$(echo "scale=1; 100 - ($OUTPUT_BYTES * 100 / $INPUT_BYTES)" | bc)

echo -e "${GREEN}→${NC} Size reduction: ${BOLD}${REDUCTION}%${NC}"
echo ""
echo -e "${YELLOW}NEXT STEPS:${NC}"
echo "  1. Update the sampler section in config.yml to use the new file:"
echo -e "     ${CYAN}input_file: data/input/lombardy.osh.pbf${NC}"
echo ""
echo "  2. Run the sampler:"
echo -e "     ${CYAN}bash run_sampler.sh${NC}"
echo ""
echo "================================================================================"
