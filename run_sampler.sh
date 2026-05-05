#!/bin/bash
# OSM Sampler Runner
# Creates/activates conda environment and runs the sampler

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Unicode symbols
CHECK="${GREEN}✓${NC}"
CROSS="${RED}✗${NC}"
ARROW="${CYAN}→${NC}"
INFO="${BLUE}ℹ${NC}"
WARN="${YELLOW}⚠${NC}"

echo ""
echo "================================================================================"
echo -e "${BOLD}${MAGENTA}╔═══════════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}${MAGENTA}║${NC}              ${BOLD}OSM ESTABLISHMENT SAMPLER - SETUP & RUN${NC}                  ${BOLD}${MAGENTA}║${NC}"
echo -e "${BOLD}${MAGENTA}╚═══════════════════════════════════════════════════════════════════════════╝${NC}"
echo "================================================================================"
echo ""

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Configuration
ENV_NAME="osm_sampler_env"
CONFIG_FILE="config.yml"

echo -e "${BOLD}[STEP 1/4] Checking Dependencies${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

# Check if conda is available
echo -ne "  Checking for conda... "
if ! command -v conda &> /dev/null; then
    echo -e "${CROSS}"
    echo -e "${RED}${BOLD}ERROR:${NC} conda is not installed or not in PATH"
    echo -e "${INFO} Install from: ${CYAN}https://docs.conda.io/en/latest/miniconda.html${NC}"
    exit 1
fi
echo -e "${CHECK}"

# Initialize conda for bash
eval "$(conda shell.bash hook)"

echo ""
echo -e "${BOLD}[STEP 2/4] Setting Up Environment${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

# Check if environment exists
echo -ne "  Checking environment '${CYAN}${ENV_NAME}${NC}'... "
if conda env list | grep -q "^${ENV_NAME} "; then
    echo -e "${CHECK}"
    echo -e "  ${ARROW} Activating existing environment..."
    conda activate "${ENV_NAME}"
else
    echo -e "${WARN} Not found"
    echo -e "  ${ARROW} Creating new environment with Python 3.11..."
    echo ""
    
    # Create environment with Python 3.11
    conda create -n "${ENV_NAME}" python=3.11 -y -q
    conda activate "${ENV_NAME}"
    
    echo ""
    echo -e "  ${ARROW} Installing conda packages..."
    conda install -n "${ENV_NAME}" -c conda-forge \
        geopandas \
        pandas \
        pyyaml \
        shapely \
        pyproj \
        boost-cpp \
        -y -q
    
    echo -e "  ${ARROW} Installing pip packages..."
    pip install -q osmium tqdm
    
    echo -e "  ${CHECK} Environment created successfully"
fi

echo ""
echo -ne "  Verifying packages... "
python -c "import osmium; import geopandas; import yaml; import shapely; import pyproj" 2>/dev/null || {
    echo -e "${CROSS}"
    echo -e "${RED}${BOLD}ERROR:${NC} Package verification failed"
    exit 1
}
echo -e "${CHECK}"

echo ""
echo -e "${BOLD}[STEP 3/4] Validating Configuration${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

# Check if config file exists
echo -ne "  Config file (${CYAN}${CONFIG_FILE}${NC})... "
if [ ! -f "${CONFIG_FILE}" ]; then
    echo -e "${CROSS}"
    echo -e "${RED}${BOLD}ERROR:${NC} Config file not found"
    exit 1
fi
echo -e "${CHECK}"

# Check if input file exists
INPUT_FILE=$(python -c "import yaml; cfg=yaml.safe_load(open('${CONFIG_FILE}')); s=cfg.get('sampler', cfg); print(s.get('input_file', ''))")
echo -ne "  Input file (${CYAN}$(basename ${INPUT_FILE})${NC})... "
if [ ! -f "${INPUT_FILE}" ]; then
    echo -e "${CROSS}"
    echo ""
    echo -e "  ${WARN} ${YELLOW}Input file not found:${NC} ${INPUT_FILE}"
    echo -e "  ${INFO} Place your OSM history file (.osh.pbf) in: ${CYAN}data/input/${NC}"
    echo -e "  ${INFO} Or update 'input_file' in the sampler section of ${CONFIG_FILE}"
    echo ""
    read -p "  Continue anyway? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${CROSS} Aborted"
        exit 1
    fi
else
    echo -e "${CHECK}"
    # Get file size
    FILE_SIZE=$(du -h "${INPUT_FILE}" | cut -f1)
    echo -e "  ${ARROW} File size: ${CYAN}${FILE_SIZE}${NC}"
fi

# Display config info
echo ""
echo -e "  ${BOLD}Configuration Summary:${NC}"
CENTER_LAT=$(python -c "import yaml; cfg=yaml.safe_load(open('${CONFIG_FILE}')); s=cfg.get('sampler', cfg); print(s.get('center_lat', ''))")
CENTER_LON=$(python -c "import yaml; cfg=yaml.safe_load(open('${CONFIG_FILE}')); s=cfg.get('sampler', cfg); print(s.get('center_lon', ''))")
RADIUS=$(python -c "import yaml; cfg=yaml.safe_load(open('${CONFIG_FILE}')); s=cfg.get('sampler', cfg); print(s.get('radius_meters', ''))")
START_DATE=$(python -c "import yaml; cfg=yaml.safe_load(open('${CONFIG_FILE}')); s=cfg.get('sampler', cfg); print(s.get('start_date', ''))")
END_DATE=$(python -c "import yaml; cfg=yaml.safe_load(open('${CONFIG_FILE}')); s=cfg.get('sampler', cfg); print(s.get('end_date', ''))")

echo -e "  ${ARROW} Center: ${CYAN}${CENTER_LAT}, ${CENTER_LON}${NC}"
echo -e "  ${ARROW} Radius: ${CYAN}${RADIUS}m${NC} ($(python -c "print(f'{${RADIUS}/1000:.1f}')") km)"
echo -e "  ${ARROW} Date range: ${CYAN}${START_DATE}${NC} to ${CYAN}${END_DATE}${NC}"

echo ""
echo -e "${BOLD}[STEP 4/4] Running Sampler${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

# Run the sampler
python scripts/osm_sampler.py --config "${CONFIG_FILE}"

SAMPLER_EXIT_CODE=$?

if [ $SAMPLER_EXIT_CODE -eq 0 ]; then
    echo ""
    echo -e "${GREEN}${BOLD}╔═══════════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}${BOLD}║${NC}                     ${CHECK} ${BOLD}SAMPLING COMPLETED SUCCESSFULLY${NC}                    ${GREEN}${BOLD}║${NC}"
    echo -e "${GREEN}${BOLD}╚═══════════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    echo -e "${INFO} ${BOLD}Next steps:${NC}"
    echo -e "  ${ARROW} Check output in: ${CYAN}data/output/samples/${NC}"
    echo -e "  ${ARROW} Open GeoPackage in QGIS for visualization"
    echo -e "  ${ARROW} Analyze CSV with pandas/R"
    echo ""
    echo -e "${INFO} ${BOLD}To run again:${NC}"
    echo -e "  ${CYAN}bash run_sampler.sh${NC}"
    echo ""
    echo -e "${INFO} ${BOLD}To activate environment manually:${NC}"
    echo -e "  ${CYAN}conda activate ${ENV_NAME}${NC}"
    echo ""
else
    echo ""
    echo -e "${RED}${BOLD}╔═══════════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}${BOLD}║${NC}                       ${CROSS} ${BOLD}SAMPLING FAILED${NC}                              ${RED}${BOLD}║${NC}"
    echo -e "${RED}${BOLD}╚═══════════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    echo -e "${INFO} Check the error messages above for details"
    echo ""
    exit 1
fi
