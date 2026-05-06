# Thesis Codebase: OSM-Based Supermarket Impact Analysis Pipeline

A comprehensive R/Python pipeline for analyzing the geographic and temporal impact of supermarket openings using OpenStreetMap data, difference-in-differences (DiD) methodology with stacked event-study design, and Poisson Panel Maximum Likelihood (PPML) estimation.

## Overview

This pipeline processes OSM data for a geographic region (Lombardy), creates spatial treatments (Area of Treatment around supermarket centroids), generates time-series panels of establishment counts, and performs robust causal inference analysis using modern DiD estimators with heterogeneity-robust inference and sensitivity analyses.

**Key Features:**
- OSM data extraction and cleaning for Italian regions
- Spatial buffering and establishment-POI/POC geocoding
- Long-format panel construction with treatment cohorts
- Stacked PPML event-study estimation with fixed effects
- HonestDiD sensitivity bounds under alternative parallel trends assumptions
- Pre-trend testing and anticipation robustness checks
- Comprehensive metadata tracking and reproducible results

## Project Structure

```
EXPORT/
├── README.md                              # This file
├── config.yml                             # Master configuration file
│
├── scripts/                               # Core analysis scripts
│   ├── osm_sampler.py                     # Extract/filter OSM data by region
│   ├── panel_builder.py                   # Generate treatment AOTs and panels
│   ├── stacked_ppml_event_study.R         # Main PPML analysis orchestrator
│   └── install_honestdid.R                # HonestDiD dependency installer
│
├── run_sampler.sh                         # Wrapper: OSM extraction (macOS/Linux)
├── run_sampler.bat                        # Wrapper: OSM extraction (Windows)
├── run_panel_builder.sh                   # Wrapper: Panel generation (macOS/Linux)
├── run_panel_builder.bat                  # Wrapper: Panel generation (Windows)
├── run_stacked_ppml_event_study.sh        # Wrapper: Main PPML analysis (macOS/Linux)
├── run_stacked_ppml_event_study.bat       # Wrapper: Main PPML analysis (Windows)
├── extract_lombardy.sh                    # Quick Lombardy extraction
│
└── data/
    ├── input/
    │   ├── .gitkeep
│   │   └── lombardy.osh.pbf               # Downloaded automatically by the sampler if missing
    └── output/
        ├── .gitkeep
        ├── samples/                       # OSM sampled data (sampler outputs)
        │   └── .gitkeep
        ├── panels/                        # Generated treatment panels
        │   └── .gitkeep
        └── analysis/                      # Analysis results
            └── .gitkeep


```

## Configuration

All workflows are controlled via `config.yml`, which defines:

### Sampler Section
```yaml
sampler:
  aoi_center:
    lon: 9.1895
    lat: 45.4642
  aoi_radius_km: 100
  data_source: "osmium"
```

### Panel Section
```yaml
panel:
  aot_radius_m: 500              # Buffer radius around supermarket centroids
  time_aggregation: quarterly     # Aggregation unit (quarterly, monthly, yearly)
  treatment_type: "supermarket"   # POI category
  outcome_types:                  # Establishment categories to count
    - "retail"
    - "food"
```

### Analysis Section
```yaml
analysis:
  pre_periods: 12                 # Quarters before treatment
  post_periods: 12                # Quarters after treatment
  anticipation:
    enabled: true                 # Run shift analysis for robustness
    periods: 6                    # Quarters to shift treatment forward
  clustering: "site_id"           # Clustering variable for SE
  fe_variables:
    - "site_stack"
    - "time_stack"
```

## Quick Start

### 1. Extract OSM Data for Your Region

```bash
# macOS/Linux:
# Edit config.yml to set your region (aoi_center, aoi_radius_km)
./run_sampler.sh
# Output: data/output/samples/<timestamp>/
```

**Windows:**
```batch
REM Edit config.yml to set your region (aoi_center, aoi_radius_km)
run_sampler.bat
REM Output: data\output\samples\<timestamp>\
```

If `data/input/lombardy.osh.pbf` is missing, the sampler downloads it automatically from Google Drive using `sampler.google_drive_file_id` in `config.yml` or the `OSM_GOOGLE_DRIVE_FILE_ID` environment variable.

### 2. Generate Treatment Panels

```bash
# macOS/Linux:
# Uses latest sample and config AOT radius
./run_panel_builder.sh
# Output: data/output/panels/panel_<timestamp>/panel.csv
```

**Windows:**
```batch
REM Uses latest sample and config AOT radius
run_panel_builder.bat
REM Output: data\output\panels\panel_<timestamp>\panel.csv
```

### 3. Run PPML Event-Study Analysis

```bash
# macOS/Linux:
# Uses latest panel, estimates stacked PPML models
./run_stacked_ppml_event_study.sh
# Output: data/output/analysis/stacked_ppml_es_<timestamp>_baseline_<year_range>/
```

**Windows:**
```batch
REM Uses latest panel, estimates stacked PPML models
run_stacked_ppml_event_study.bat
REM Output: data\output\analysis\stacked_ppml_es_<timestamp>_baseline_<year_range>\
```

### 4. Run Complete Pipeline

```bash
# macOS/Linux:
# Execute sampler → panel builder → analysis in sequence
./run_complete_analysis.sh
```

## Detailed Workflow

### Step 1: OSM Sampling (`osm_sampler.py`)

**Input:**
- OSM PBF file (e.g., `data/input/lombardy.osh.pbf`) downloaded from Google Drive, originally sourced from Geofabrik
- AOI center (lon, lat) and radius (km) from `config.yml`

**Process:**
- Extracts establishments within circular AOI using osmium CLI
- Filters by OSM tags (e.g., `shop=supermarket` for treatment)
- Standardizes columns and removes duplicates
- Outputs GeoDataFrame with geometry

**Output:**
- CSV with columns: `osm_id, name, lat, lon, shop, geometry, extracted_at`

**Key Parameters:**
- `aoi_center.lon`, `aoi_center.lat`: AOI center coordinates
- `aoi_radius_km`: Radius of analysis region

---

### Step 2: Panel Generation (`panel_builder.py`)

**Input:**
- Establishments from sampler output
- AOT radius (m) from `config.yml`
- Time range for analysis

**Process:**
1. **Create Area of Treatment (AOT):**
   - Buffer each treatment unit (supermarket) by `aot_radius_m`
   - Create circular AOTs as GeoDataFrame

2. **Spatial Join:**
   - Intersect outcomes (e.g., retail, food establishments) with AOTs
   - Each outcome linked to nearest treatment unit

3. **Time-Series Panel:**
   - Aggregate establishment counts by site, month (or quarter/year)
   - Calculate time-to-treatment for each site
   - Assign cohort based on treatment timing
   - Filter to balanced window (pre/post periods from config)

4. **Output Format (Long):**
   ```
   site_id, month, poi_poc, count, entries, exits, 
   date_of_treatment, cohort_id, time_to_treatment
   ```

**Output:**
- `data/output/panels/panel_<timestamp>/panel.csv`
- Panel statistics printed to console

**Key Parameters:**
- `aot_radius_m`: Treatment buffer radius
- `time_aggregation`: "monthly", "quarterly", or "yearly"
- `treatment_type`: POI category (e.g., "supermarket")
- `outcome_types`: Array of establishment types to count

---

### Step 3: PPML Event-Study Analysis (`stacked_ppml_event_study.R`)

**Input:**
- Latest panel from `data/output/panels/`
- Config settings for pre/post periods, clustering, anticipation

**Process:**

1. **Read and Validate Panel:**
   - Load latest panel automatically
   - Print metadata: sampler AOI, panel AOT radius, treatment type, outcomes

2. **Construct Stacked Dataset:**
   - Create separate cohort stacks for each treatment timing
   - Include not-yet-treated units as controls
   - Follows Callaway & Sant'Anna (2021)

3. **Estimate PPML Models:**
   - **Baseline Model:**
     ```
     fepois(count ~ i(time_to_treatment, ref=-1) | site_stack + time_stack, cluster="site_id")
     ```
   - Event-study coefficients give Average Treatment Effect on Treated (ATT) per period
   - Fixed effects absorb site-level and time-level heterogeneity

4. **Pre-Trend Testing:**
   - Wald test on all pre-treatment coefficients (H₀: β_pre = 0)
   - Printed to results file

5. **HonestDiD Sensitivity:**
   - Computes bounds on ATT under alternative parallel trends assumptions
   - Tests robustness to user-specified violations of parallel trends
   - Outputs confidence sets (eff-balanced, H-optimized, and ARP)

6. **Anticipation Analysis (Optional):**
   - If `anticipation.enabled: true`:
     - Shifts treatment indicator forward by `anticipation.periods`
     - Re-estimates models to test for pre-treatment effects
     - Compares coefficients to detect anticipatory behavior

**Output:**
- `data/output/analysis/stacked_ppml_es_<timestamp>_baseline_<year_range>/results.txt`
  - Model summaries, pre-trend test, coefficient estimates
- `.../_shift_<n>quarters/results.txt` (if anticipation enabled)
  - Shifted analysis results
- `.../_by_outcome_type/` (optional)
  - Separate estimates by establishment type
- PNG plots of event-study coefficients

**Key Parameters:**
- `pre_periods`: Quarters before treatment to include
- `post_periods`: Quarters after treatment to include
- `clustering`: Variable for robust SE (usually "site_id")
- `fe_variables`: Fixed effect specifications
- `anticipation.enabled`: Enable shift analysis
- `anticipation.periods`: Number of periods to shift forward

---

## Configuration Guide

### `config.yml` Master Configuration

All scripts read from `config.yml`. Update this single file to control the entire pipeline.

#### Sampler Section
```yaml
sampler:
  aoi_center:
    lon: 9.1895      # Longitude of AOI center
    lat: 45.4642     # Latitude of AOI center
  aoi_radius_km: 100 # Radius of analysis region (km)
  data_source: "osmium"
  google_drive_file_id: "" # Optional: centralized download ID for the OSM file
```

#### Panel Section
```yaml
panel:
  aot_radius_m: 500                    # Treatment buffer radius (meters)
  time_aggregation: quarterly           # Time unit: quarterly, monthly, yearly
  treatment_type: "supermarket"         # Primary treatment unit OSM tag
  outcome_types:                        # Establishment types to count
    - "retail"
    - "food"
```

#### Analysis Section
```yaml
analysis:
  pre_periods: 12                       # Pre-treatment quarters
  post_periods: 12                      # Post-treatment quarters
  clustering: "site_id"                 # Variable for clustered SE
  fe_variables:
    - "site_stack"                      # Site fixed effects
    - "time_stack"                      # Time fixed effects
  anticipation:
    enabled: true                       # Run anticipation checks
    periods: 6                          # Shift quarters forward
```

### Adjusting Analysis Parameters

- **Change region:** Update `aoi_center` and `aoi_radius_km` in `config.yml`, then run `./run_sampler.sh`
- **Use a different OSM file:** Update `sampler.google_drive_file_id` or set `OSM_GOOGLE_DRIVE_FILE_ID`
- **Change treatment buffer:** Update `aot_radius_m` in `config.yml`, then run `./run_panel_builder.sh`
- **Change time window:** Update `pre_periods` and `post_periods`, then re-run analysis
- **Check anticipation:** Set `anticipation.enabled: true` and adjust `anticipation.periods` (default 6 quarters)

## Running Individual Steps

### Extract OSM Data Only
```bash
./run_sampler.sh
```

### Generate Panel Only
```bash
# Prerequisites: Latest sample in data/output/samples/
./run_panel_builder.sh
```

### Run Analysis Only
```bash
# Prerequisites: Latest panel in data/output/panels/
./run_stacked_ppml_event_study.sh
```

### Quick Lombardy Extraction
```bash
# Shortcut for extracting Lombardy OSM data with preset coordinates
./extract_lombardy.sh
```

## Output Files

### Analysis Results Directory Structure
```
data/output/analysis/stacked_ppml_es_<timestamp>_baseline_<year_range>/
├── results.txt                    # Main results (model summaries, coefficients)
├── event_study_plot.png           # Event-study coefficient plot
├── honestdid_bounds.png           # HonestDiD sensitivity bounds plot
├── shift_<n>quarters/
│   ├── results.txt                # Shifted analysis results
│   └── event_study_plot.png
└── by_outcome_type/
    ├── <outcome_type>_results.txt
    └── <outcome_type>_plot.png
```

### Key Output Files

- **results.txt**: Contains:
  - Sampler metadata (AOI, radius)
  - Panel metadata (AOT radius, time range, cohort distribution)
  - PPML model summary (coefficients, SE, p-values)
  - Pre-trend test results (Wald statistic, p-value)
  - HonestDiD bounds (lower/upper bounds on ATT)
  - Interpretation guidance

- **event_study_plot.png**: Coefficient plot with 95% CI across event-study window

- **honestdid_bounds.png**: Sensitivity bounds under parallel trends violations

## Dependencies

### R Packages
- `fixest` (PPML estimation)
- `tidyverse` (dplyr, ggplot2, tidyr)
- `data.table` (fast data manipulation)
- `sf` (spatial operations)
- `honestdid` (sensitivity analysis)
- `yaml` (config reading)

### Python Packages
- `geopandas` (spatial data)
- `pandas` (data manipulation)
- `pyyaml` (YAML config)
- `shapely` (geometry operations)

### System Requirements
- R >= 4.0
- Python >= 3.8
- osmium (CLI tool for OSM data extraction)
- git (for cloning)

### Installation

```bash
# R dependencies (from R console)
install.packages(c("fixest", "tidyverse", "data.table", "sf", "yaml"))

# HonestDiD (special installation)
source("scripts/install_honestdid.R")

# Python dependencies
pip install geopandas pandas pyyaml shapely

# System tools (macOS)
brew install osmium-tool
```

## Windows Setup & Usage Guide

The pipeline includes native Windows batch file wrappers (`.bat` files) for all three main steps, providing the same functionality and experience as the macOS/Linux shell scripts.

### Windows Prerequisites

**1. Install Python 3.10+**
- Download from https://www.python.org/downloads/
- **Important:** Check "Add Python to PATH" during installation
- Verify: Open Command Prompt and run `python --version`

**2. Install R 4.0+**
- Download from https://cran.r-project.org/
- Accept default installation options
- Verify: Open Command Prompt and run `Rscript --version`

**3. Install Required Python Packages**
```batch
pip install geopandas pandas pyyaml shapely osmium tqdm
```

**4. Install OSM Extraction Tool (Osmium)**
```batch
pip install osmium
```

**5. Install R Packages**
Open Command Prompt and run:
```batch
Rscript -e "install.packages(c('fixest', 'ggplot2', 'dplyr', 'tidyr', 'data.table', 'yaml', 'zoo', 'car'))"
```

Then install HonestDiD:
```batch
Rscript -e "if (!require('devtools')) install.packages('devtools'); devtools::install_github('asheshrambachan/HonestDiD')"
```

### Running the Pipeline on Windows

Navigate to the EXPORT folder in Command Prompt and run the batch files:

**1. Extract OSM Data**
```batch
run_sampler.bat
```
- Validates config.yml
- Checks Python dependencies
- Downloads OSM file from Google Drive if missing
- Extracts establishments for your region
- Output: `data\output\samples\<timestamp>\`

**2. Generate Treatment Panels**
```batch
run_panel_builder.bat
```
- Auto-detects latest sample
- Creates Area of Treatment (AOT) buffers
- Generates time-series panel
- Output: `data\output\panels\panel_<timestamp>\panel.csv`

Optionally specify a sample directory:
```batch
run_panel_builder.bat data\output\samples\sample_name
```

**3. Run PPML Event-Study Analysis**
```batch
run_stacked_ppml_event_study.bat
```
- Checks R installation and packages
- Auto-installs missing dependencies
- Validates panel data
- Runs stacked PPML estimation with HonestDiD
- Output: `data\output\analysis\stacked_ppml_es_<timestamp>_baseline_<year_range>\`

### Batch File Features

All Windows batch files include:

- **Dependency checking:** Validates Python, R, and required packages
- **Auto-installation:** Missing R packages are installed automatically
- **Error handling:** Clear error messages if prerequisites are missing
- **Progress tracking:** Visual feedback with [OK], [FAIL], and [>>] indicators
- **File validation:** Checks for required input files before running
- **Pause on exit:** Results remain visible after script completion

### Windows-Specific Notes

1. **File Paths:** Use backslashes (`\`) in batch files, but Python/R handle both `/` and `\`
2. **Command Prompt:** Run batch files from Command Prompt or PowerShell
3. **Conda Users:** If using Conda, batch files use the system Python by default. To use a Conda environment:
   ```batch
   conda activate your_env_name
   run_sampler.bat
   ```
4. **Spaces in Paths:** All batch files handle paths with spaces correctly
5. **Output Directories:** Windows batch files create the same output structure as Unix scripts

### Troubleshooting Windows Issues

**"Python is not installed or not in PATH"**
- Reinstall Python and ensure "Add Python to PATH" is checked
- Restart Command Prompt after installation
- Verify: `python --version`

**"R is not installed or not in PATH"**
- Reinstall R with default options
- Verify: `Rscript --version`

**"Package installation failed"**
- Run batch files as Administrator (right-click → Run as Administrator)
- Or manually install packages:
  ```batch
  Rscript -e "install.packages('package_name', repos = 'https://cloud.r-project.org/')"
  ```

**"No panel directories found"**
- Ensure `run_sampler.bat` completed successfully
- Check that `data\output\samples\` contains files
- Check for error messages in the sampler output

**Permission Denied Errors**
- Right-click the Command Prompt and select "Run as Administrator"
- Alternatively, use PowerShell with admin privileges

## Troubleshooting

### "No panel directories found"
- Ensure `./run_panel_builder.sh` has been run successfully
- Check that `data/output/panels/` contains `panel_*/panel.csv` files

### "Could not find latest sample"
- Ensure `./run_sampler.sh` has been run successfully
- Check that `data/output/samples/` contains output directories

### PPML convergence issues
- Check panel for large count values; may need aggregation
- Verify treatment variation (some sites need treatment exposure)
- Inspect pre-trend test; if significant, consider different specification

### Spatial join producing too few observations
- Check AOT radius is reasonable (500m = default)
- Verify treatment and outcome OSM tags are present in data
- Run diagnostic plots from panel builder

## Citation & References

**Core Methodology:**
- Callaway, B., & Sant'Anna, P. H. (2021). "Difference-in-Differences with Multiple Time Periods." *Journal of Econometrics*, 225(2), 200-230.

**Sensitivity Analysis:**
- Rambachan, A., & Roth, J. (2023). "Honest Inference on Causal Effects Under Identification Uncertainty." *Econometric Reviews*, 42(3-4), 226-265.

**OSM Data:**
- OpenStreetMap contributors. "OpenStreetMap." https://www.openstreetmap.org/
- Geofabrik. "Download OpenStreetMap Data." https://www.geofabrik.de/

## Contact & Support

For questions, issues, or feature requests, please file an issue or contact the author.

---

**Last Updated:** May 2026  
**Status:** Production-Ready
