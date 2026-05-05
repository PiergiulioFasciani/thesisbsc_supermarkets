#!/bin/bash

# ==============================================================================
# Run Stacked PPML Event-Study Analysis with Dependency Checks
# ==============================================================================
# This script:
# 1. Checks for required R packages and installs missing ones
# 2. Verifies HonestDiD is properly installed
# 3. Runs the stacked PPML event-study analysis
# 4. Reports the results location
#
# Usage:
#   ./run_stacked_ppml_event_study.sh
# ======================================================================

set -e  # Exit on error

echo "═══════════════════════════════════════════════════════════════"
echo "  STACKED PPML EVENT-STUDY ANALYSIS - DEPENDENCY CHECK & EXECUTION"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# ==============================================================================
# 1. Check R Installation
# ==============================================================================

echo "─────────────────────────────────────────────────────────────"
echo "1. Checking R installation..."
echo "─────────────────────────────────────────────────────────────"

if ! command -v Rscript &> /dev/null; then
    echo "ERROR: R is not installed or not in PATH"
    echo "Please install R from: https://cran.r-project.org/"
    exit 1
fi

R_VERSION=$(Rscript --version 2>&1 | head -n 1)
echo "✓ R found: $R_VERSION"
echo ""

# ==============================================================================
# 2. Check Required R Packages
# ==============================================================================

echo "─────────────────────────────────────────────────────────────"
echo "2. Checking required R packages..."
echo "─────────────────────────────────────────────────────────────"

REQUIRED_PACKAGES=(
    "fixest"
    "ggplot2"
    "dplyr"
    "tidyr"
    "data.table"
    "yaml"
    "zoo"
    "car"
    "HonestDiD"
)

MISSING_PACKAGES=()

for package in "${REQUIRED_PACKAGES[@]}"; do
    if Rscript -e "if (!require('$package', quietly = TRUE)) quit(status = 1)" &> /dev/null; then
        echo "  ✓ $package"
    else
        echo "  ✗ $package (missing)"
        MISSING_PACKAGES+=("$package")
    fi
done

echo ""

# ==============================================================================
# 3. Install Missing Packages
# ==============================================================================

if [ ${#MISSING_PACKAGES[@]} -gt 0 ]; then
    echo "─────────────────────────────────────────────────────────────"
    echo "3. Installing missing packages..."
    echo "─────────────────────────────────────────────────────────────"

    for package in "${MISSING_PACKAGES[@]}"; do
        echo "Installing $package..."

        if [ "$package" = "HonestDiD" ]; then
            echo "  Installing HonestDiD from GitHub (requires devtools)..."
            Rscript -e "
                if (!require('devtools', quietly = TRUE)) {
                    install.packages('devtools', repos = 'https://cloud.r-project.org/')
                }
                devtools::install_github('asheshrambachan/HonestDiD')
            "
        else
            Rscript -e "install.packages('$package', repos = 'https://cloud.r-project.org/')"
        fi

        if [ $? -eq 0 ]; then
            echo "  ✓ $package installed successfully"
        else
            echo "  ✗ ERROR: Failed to install $package"
            exit 1
        fi
    done

    echo ""
    echo "All missing packages installed successfully!"
    echo ""
else
    echo "All required packages are already installed!"
    echo ""
fi

# ==============================================================================
# 4. Verify Configuration File
# ==============================================================================

echo "─────────────────────────────────────────────────────────────"
echo "4. Verifying configuration..."
echo "─────────────────────────────────────────────────────────────"

CONFIG_FILE="config.yml"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Configuration file not found: $CONFIG_FILE"
    echo "Please create the configuration file before running analysis"
    exit 1
fi

echo "✓ Configuration file found: $CONFIG_FILE"

echo ""
echo "Configuration settings:"
echo "─────────────────────────────────────────────────────────────"
grep -E "(time_aggregation|estimator|start_date|end_date|pre_periods|post_periods|cluster_var)" "$CONFIG_FILE" | sed 's/^/  /'
echo "─────────────────────────────────────────────────────────────"
echo ""

# ==============================================================================
# 5. Check for Panel Data
# ==============================================================================

echo "─────────────────────────────────────────────────────────────"
echo "5. Checking for panel data..."
echo "─────────────────────────────────────────────────────────────"

PANEL_DIR="data/output/panels"

if [ ! -d "$PANEL_DIR" ]; then
    echo "ERROR: Panel directory not found: $PANEL_DIR"
    echo "Please run panel_builder.py first to create panel data"
    exit 1
fi

PANEL_COUNT=$(find "$PANEL_DIR" -name "panel.csv" | wc -l | tr -d ' ')

if [ "$PANEL_COUNT" -eq 0 ]; then
    echo "ERROR: No panel.csv files found in $PANEL_DIR"
    echo "Please run panel_builder.py first to create panel data"
    exit 1
fi

LATEST_PANEL=$(python3 - <<'PY'
from pathlib import Path
panel_files = [p for p in Path('data/output/panels').glob('panel_*/panel.csv') if p.is_file()]
if panel_files:
    latest = max(panel_files, key=lambda p: p.stat().st_mtime)
    print(latest)
PY
)
echo "✓ Found $PANEL_COUNT panel file(s)"
echo "  Latest panel: $LATEST_PANEL"
echo ""

# ==============================================================================
# 6. Run Stacked PPML Event-Study Analysis
# ==============================================================================

echo "═══════════════════════════════════════════════════════════════"
echo "  RUNNING STACKED PPML EVENT-STUDY ANALYSIS"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "This will perform:"
echo "  - Standard stacked PPML event studies for count, entries, exits"
echo "  - Wald pretrend tests on the closest pre-treatment periods"
echo "  - Support-weighted post-treatment PPML ATT summaries"
echo "  - HonestDiD sensitivity analysis for all specifications"
echo ""
echo "Starting analysis..."
echo ""

START_TIME=$(date +%s)

Rscript scripts/stacked_ppml_event_study.R

EXIT_CODE=$?
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$((ELAPSED / 60))
ELAPSED_SEC=$((ELAPSED % 60))

echo ""
echo "═══════════════════════════════════════════════════════════════"

if [ $EXIT_CODE -eq 0 ]; then
    echo "  ✓ ANALYSIS COMPLETED SUCCESSFULLY"
    echo "═══════════════════════════════════════════════════════════════"
    echo ""
    echo "Execution time: ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
    echo ""

    OUTPUT_DIR=$(find data/output/analysis -type d -name "stacked_ppml_es_*" | sort -r | head -n 1)

    if [ -n "$OUTPUT_DIR" ]; then
        echo "Results location:"
        echo "  Directory: $OUTPUT_DIR"
        echo ""
        echo "Output files:"
        echo "  - stacked_ppml_event_study_results.txt (Full results text file)"
        echo "  - csv/complete_summary.csv             (Summary table)"
        echo "  - plots/standard/                     (Event-study and HonestDiD PNGs)"
        echo ""

        PLOT_COUNT=$(find "$OUTPUT_DIR/plots" -name "*.png" | wc -l | tr -d ' ')
        echo "  Generated $PLOT_COUNT plots"
        echo ""

        SUMMARY_FILE="$OUTPUT_DIR/csv/complete_summary.csv"
        if [ -f "$SUMMARY_FILE" ]; then
            echo "Quick summary (see full results in stacked_ppml_event_study_results.txt):"
            echo "─────────────────────────────────────────────────────────────"
            head -n 10 "$SUMMARY_FILE"
            echo "─────────────────────────────────────────────────────────────"
        fi
    fi
else
    echo "  ✗ ANALYSIS FAILED"
    echo "═══════════════════════════════════════════════════════════════"
    echo ""
    echo "Exit code: $EXIT_CODE"
    echo "Execution time: ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
    echo ""
    echo "Please check the error messages above for details"
    exit $EXIT_CODE
fi

echo ""
echo "═══════════════════════════════════════════════════════════════"
echo "  COMPLETE"
echo "═══════════════════════════════════════════════════════════════"
