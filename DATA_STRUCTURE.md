# DDD Analysis Data Structure

This document describes the data organization for the Difference-in-Difference-in-Differences (DDD) analysis.

## Directory Structure

```
data/
├── pbf/                      # OSM PBF input files (original data)
├── raw/                      # Raw DDD extraction outputs
│   └── ddd_sample_TIMESTAMP_RADIUS/
│       ├── ddd_sample.csv    # Feature-level data with timestamps
│       └── ddd_sample.gpkg   # GeoPackage with spatial layers
└── processed/                # Processed panel data
    └── ddd_panel_TIMESTAMP_RADIUS.csv  # Panel dataset for regression
```

## Pipeline Workflow

### Step 1: Extract Raw Data (`main_ddd.py`)

**Input:** OSM PBF file + center coordinates + radius  
**Output:** `data/raw/ddd_sample_*/`

Extracts:
- **Supermarkets** within AOI + 1000m buffer
- **Treatment POIs** (bakery, butcher, greengrocer, convenience) in AOI + buffer
- **Control POCs** (funeral_home, hairdresser, restaurant) in AOI only
- All features include OSM timestamps and distance to Duomo

**GPKG Layers:**
1. `pois` - Treatment points of interest
2. `pocs` - Points of control (negative controls)
3. `supermarkets` - Supermarket locations
4. `duomo` - Reference point (marked as star)
5. `aoi_circle` - Area of Interest boundary
6. `buffer_circle` - Buffer zone boundary

### Step 2: Process into Panel (`process_ddd_panel.py`)

**Input:** `data/raw/ddd_sample_*/`  
**Output:** `data/processed/ddd_panel_*.csv`

Creates a balanced panel dataset with:
- Site-level clustering of supermarkets (75m radius)
- Ring-based spatial treatment (in: 0-500m, mid: 500-1000m, out: 1000-1500m)
- Monthly time series (2018-01 to 2025-09, 93 months)
- Stock-based outcome counts
- Treatment and exposure indicators

## Panel Data Structure

### Primary Keys
- **area_id**: Unique micro-area (format: `SXXXX_in/mid/out`)
- **site_id**: Supermarket site identifier (format: `SXXXX`)
- **month**: Monthly time index (YYYY-MM-01)
- **category**: Business category (canonical names)

### Treatment Variables
- **A** (0/1): Affected indicator (1 = treatment, 0 = placebo)
- **t0** (date): Site entry month (earliest supermarket opening)
- **E** (0/1): Exposure indicator (1 if month ≥ t0)
- **event_time** (int): Months since entry (negative pre, positive post)
- **ev_cap** (int): Capped event time [-24, +24]

### Spatial Design
- **ring** (str): Distance band (in/mid/out)
- **ring_min_m**, **ring_max_m** (int): Ring boundaries in meters
- **area_km2** (float): Ring area in km²

### Outcome Variable
- **y** (int ≥ 0): Stock count of active units (cumulative entries)

### Eligibility
- **eligible** (0/1): Site has ≥12 months pre & post within panel window

### Metadata
- **sample_window_start**, **sample_window_end** (date): Panel boundaries

## Data Integrity Rules

✅ **Uniqueness:** Each (area_id, category, month) appears exactly once  
✅ **Completeness:** All areas have all 93 months (balanced panel)  
✅ **Zero-filling:** Missing categories get y=0 (not NA)  
✅ **Constancy:** A is constant within category; t0 within site_id; ring attributes within area_id  
✅ **No NAs:** All model variables are non-missing

## Treatment Assignment

### Affected Categories (A=1)
- bakery
- butcher
- greengrocer
- convenience

### Placebo Categories (A=0)
- funeral_home (funeral_directors in raw data)
- hairdresser
- restaurant

## Usage Example

```bash
# Step 1: Extract raw data
python python/main_ddd.py \
  --center "45.464200,9.191400" \
  --radius-m 1500 \
  --pbf data/pbf/nord-ovest-latest.osm.pbf \
  --out-dir data/raw

# Step 2: Process into panel
python python/process_ddd_panel.py \
  --input data/raw/ddd_sample_06102025_211634_1500 \
  --output data/processed
```

## Panel Summary Statistics

Example output:
- **Total rows:** 54,684
- **Unique areas:** 84 (71 sites × avg 1.18 rings)
- **Categories:** 7 (4 affected + 3 placebo)
- **Time periods:** 93 months (2018-01 to 2025-09)
- **Non-zero outcomes:** ~18.5%
- **Eligible observations:** ~63%

## Regression-Ready Format

The panel data is structured for difference-in-difference-in-differences estimation:

```R
# Fixed effects specification
feols(y ~ E:A + controls | area_id + month + category, 
      data = panel, 
      cluster = ~site_id)
```

Where:
- **E×A** is the triple-difference estimator
- Fixed effects absorb area, time, and category-specific trends
- Clustering at site level accounts for spatial correlation
