#!/usr/bin/env python3
"""
OSM Panel Generator
Creates a long-format site-month panel for stacked PPML event-study estimation of supermarket entry effects.
"""

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from datetime import datetime
from pathlib import Path
import logging
import argparse
import sys
import hashlib
import yaml

# ANSI color codes
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def generate_site_id(osm_id):
    """Generate alphanumeric site ID"""
    hash_val = hashlib.md5(str(osm_id).encode()).hexdigest()[:6].upper()
    return f"S{osm_id % 1000000:06d}{hash_val}"

def create_aots(supermarket_gdf, aot_radius_m=1000):
    """Create one circular Area of Treatment (AOT) per supermarket.

    AOT definition:
    A circular buffer of fixed radius r around each supermarket.
    """

    logger.info(f"Creating one AOT per supermarket with radius {aot_radius_m}m for {len(supermarket_gdf)} supermarkets...")
    
    # Project to metric CRS (UTM zone 32N for Italy/Milan)
    supermarkets_proj = supermarket_gdf.to_crs('EPSG:32632')
    
    aots_list = []
    
    for idx, site in supermarkets_proj.iterrows():
        site_id = generate_site_id(site['osm_id'])
        
        aot_geom = site.geometry.buffer(aot_radius_m)
        aot_data = {
            'site_id': site_id,
            'treatment_date': pd.to_datetime(site['opening_date']),
            'geometry': aot_geom
        }
        aots_list.append(aot_data)

    aots_gdf = gpd.GeoDataFrame(aots_list, crs='EPSG:32632')

    logger.info(f"✓ Created {len(aots_gdf)} AOTs")

    return aots_gdf

def spatial_join_establishments(aots_gdf, establishments_gdf):
    """Perform spatial join of POIs/POCs to AOTs"""
    
    logger.info(f"Spatial join: {len(establishments_gdf)} establishments to {len(aots_gdf)} AOTs...")
    
    # Project establishments to same CRS
    establishments_proj = establishments_gdf.to_crs('EPSG:32632')
    
    # Spatial join
    joined = gpd.sjoin(establishments_proj, aots_gdf, how='inner', predicate='within')
    
    # Keep relevant columns
    result = joined[[
        'osm_id', 'type', 'opening_date', 'closing_date', 'site_id', 'treatment_date'
    ]].copy()
    
    result.rename(columns={'osm_id': 'establishment_osm_id'}, inplace=True)
    
    logger.info(f"✓ Found {len(result)} establishment-AOT associations")
    
    return result

def generate_panel(spatial_join_df, start_date='2010-01-01', end_date='2025-12-31'):
    """Generate long-format panel data"""
    
    logger.info("Generating panel dataset...")
    
    # Create date range (monthly observations)
    date_range = pd.date_range(start=start_date, end=end_date, freq='MS')
    
    required_columns = [
        'site_id', 'month', 'poi_poc', 'count', 'entries', 'exits',
        'date_of_treatment', 'cohort_id', 'time_to_treatment'
    ]

    # Get unique site combinations
    sites = spatial_join_df[['site_id', 'treatment_date']].drop_duplicates()
    
    panel_rows = []
    
    for _, site_row in sites.iterrows():
        site_id = site_row['site_id']
        treatment_date = site_row['treatment_date']
        
        # Get establishments in this site
        establishments = spatial_join_df[
            (spatial_join_df['site_id'] == site_id)
        ]
        
        # For each month in the panel
        for month in date_range:
            month_str = month.strftime('%Y-%m')
            
            # Count POIs and POCs active in this month
            for poi_type in ['POI', 'POC']:
                type_establishments = establishments[establishments['type'] == poi_type]
                
                # Count active establishments (opened before/during month, not closed or closed after month)
                active = type_establishments[
                    (pd.to_datetime(type_establishments['opening_date']) <= month) &
                    ((type_establishments['closing_date'].isna()) | 
                     (pd.to_datetime(type_establishments['closing_date']) > month))
                ]
                
                count = len(active)
                
                # Count entries (opened in this month)
                entries = type_establishments[
                    pd.to_datetime(type_establishments['opening_date']).dt.to_period('M') == month.to_period('M')
                ]
                
                # Count exits (closed in this month)
                exits = type_establishments[
                    (~type_establishments['closing_date'].isna()) &
                    (pd.to_datetime(type_establishments['closing_date']).dt.to_period('M') == month.to_period('M'))
                ]
                
                # Calculate time to treatment
                months_to_treatment = (month.year - treatment_date.year) * 12 + (month.month - treatment_date.month)
                # Calculate cohort (year and quarter of treatment)
                cohort_year = treatment_date.year
                cohort_quarter = (treatment_date.month - 1) // 3 + 1
                cohort_id = f"{cohort_year}Q{cohort_quarter}"
                
                panel_row = {
                    'site_id': site_id,
                    'month': month_str,
                    'poi_poc': poi_type,
                    'count': count,
                    'entries': len(entries),
                    'exits': len(exits),
                    'date_of_treatment': treatment_date.strftime('%Y-%m'),
                    'cohort_id': cohort_id,
                    'time_to_treatment': int(months_to_treatment)
                }
                
                panel_rows.append(panel_row)
    
    panel_df = pd.DataFrame(panel_rows)
    panel_df = panel_df[required_columns]

    final_columns = panel_df.columns.tolist()
    print(f"Final panel columns: {final_columns}")

    extra_cols = sorted(set(final_columns) - set(required_columns))
    missing_cols = sorted(set(required_columns) - set(final_columns))
    if extra_cols or missing_cols:
        raise ValueError(
            f"Panel schema mismatch. Missing: {missing_cols if missing_cols else 'None'}; "
            f"Extra: {extra_cols if extra_cols else 'None'}"
        )
    
    logger.info(f"✓ Generated panel: {len(panel_df):,} observations")
    logger.info(f"  Sites: {panel_df['site_id'].nunique()}")
    logger.info(f"  Time periods: {panel_df['month'].nunique()}")
    logger.info(f"  Cohorts: {panel_df['cohort_id'].nunique()}")
    
    return panel_df

def load_config(config_path):
    """Load and parse configuration file"""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def get_workflow_config(config, section_name):
    """Return a workflow-specific section or fall back to the root config."""
    section = config.get(section_name)
    if isinstance(section, dict):
        return section
    return config

def detect_latest_sample():
    """Auto-detect most recent sample directory"""
    sample_dirs = [p for p in Path('data/output/samples').glob('sample_*') if p.is_dir()]
    sample_files = [p / 'establishments.csv' for p in sample_dirs if (p / 'establishments.csv').exists()]
    if sample_files:
        latest_sample = max(sample_files, key=lambda p: p.stat().st_mtime)
        return latest_sample
    return None

def infer_date_range(establishments_df):
    """Infer date range from establishment data"""
    all_dates = []
    
    # Collect all opening dates
    opening_dates = pd.to_datetime(establishments_df['opening_date'], errors='coerce').dropna()
    all_dates.extend(opening_dates.tolist())
    
    # Collect all closing dates
    closing_dates = pd.to_datetime(establishments_df['closing_date'], errors='coerce').dropna()
    all_dates.extend(closing_dates.tolist())
    
    if all_dates:
        min_date = min(all_dates).strftime('%Y-%m-%d')
        max_date = max(all_dates).strftime('%Y-%m-%d')
        return min_date, max_date
    
    return '2010-01-01', '2025-12-31'

def main():
    parser = argparse.ArgumentParser(
        description='Generate panel dataset from OSM establishments sample',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python scripts/panel_builder.py --config config.yml
  python scripts/panel_builder.py --input data/output/samples/sample_070226_123620/establishments.csv
        """
    )
    
    parser.add_argument('--config', '-c', default=None,
                       help='Configuration file (YAML)')
    parser.add_argument('--input', '-i',
                       help='Input establishments CSV file (overrides config)')
    parser.add_argument('--output-dir', '-o',
                       help='Output directory (overrides config)')
    parser.add_argument('--start-date',
                       help='Panel start date (overrides config)')
    parser.add_argument('--end-date',
                       help='Panel end date (overrides config)')
    
    args = parser.parse_args()
    
    direct_mode = any([args.input, args.output_dir, args.start_date, args.end_date])

    # Load configuration when requested or when running in config-driven mode.
    config = {}
    config_path = None
    if args.config:
        config_path = Path(args.config)
    elif not direct_mode:
        config_path = Path('config.yml')

    if config_path is not None:
        if not config_path.exists():
            print(f"{Colors.RED}✗ Config file not found: {config_path}{Colors.ENDC}")
            sys.exit(1)
        config = load_config(config_path)

    panel_config = get_workflow_config(config, 'panel')
    
    # Determine input file
    input_file = args.input or panel_config.get('input_file', '')
    if not input_file:
        # Auto-detect latest sample
        input_file = detect_latest_sample()
        if not input_file:
            print(f"{Colors.RED}✗ No input file specified and no samples found{Colors.ENDC}")
            sys.exit(1)
        logger.info(f"Auto-detected latest sample: {input_file}")
    
    input_file = Path(input_file)
    if not input_file.exists():
        print(f"{Colors.RED}✗ Input file not found: {input_file}{Colors.ENDC}")
        sys.exit(1)
    
    # Determine output directory
    output_dir = args.output_dir or panel_config.get('output_dir', 'data/output')
    
    # Load establishments
    print("\n" + "="*80)
    print(f"{Colors.BOLD}{Colors.HEADER}OSM PANEL GENERATOR{Colors.ENDC}")
    print("="*80)
    print(f"\n{Colors.BOLD}Configuration:{Colors.ENDC}")
    print(f"  Config file: {Colors.CYAN}{config_path if config_path else 'command-line arguments'}{Colors.ENDC}")
    print(f"  Input: {Colors.CYAN}{input_file}{Colors.ENDC}")
    print("\n" + "-"*80 + "\n")
    
    print(f"{Colors.BOLD}[1/5] Loading establishments...{Colors.ENDC}")
    establishments_df = pd.read_csv(input_file)
    logger.info(f"✓ Loaded {len(establishments_df):,} establishments")
    
    # Infer date range from data if not specified
    start_date = args.start_date or panel_config.get('start_date', '')
    end_date = args.end_date or panel_config.get('end_date', '')
    
    if not start_date or not end_date:
        inferred_start, inferred_end = infer_date_range(establishments_df)
        start_date = start_date or inferred_start
        end_date = end_date or inferred_end
        logger.info(f"✓ Inferred date range: {start_date} to {end_date}")
    
    print(f"  Date range: {Colors.CYAN}{start_date}{Colors.ENDC} to {Colors.CYAN}{end_date}{Colors.ENDC}")
    
    # Get configuration parameters
    treatment_type = panel_config.get('treatment_type', 'Supermarket')
    aot_radius_m = panel_config.get('aot_radius_m', 1000)
    outcome_types = panel_config.get('outcome_types', [])
    category_filter = panel_config.get('category_filter', {})
    min_establishments = panel_config.get('min_establishments_per_site', 0)
    treatment_start = panel_config.get('treatment_start_date', '')
    treatment_end = panel_config.get('treatment_end_date', '')
    
    # Log configuration
    logger.info(f"✓ Configuration loaded:")
    logger.info(f"  Treatment type: {treatment_type}")
    logger.info(f"  AOT radius (meters): {aot_radius_m}")
    
    # Apply filters
    if treatment_start:
        establishments_df = establishments_df[
            (establishments_df['type'] != treatment_type) |
            (pd.to_datetime(establishments_df['opening_date']) >= pd.to_datetime(treatment_start))
        ]
        logger.info(f"✓ Filtered treatments after {treatment_start}")
    
    if treatment_end:
        establishments_df = establishments_df[
            (establishments_df['type'] != treatment_type) |
            (pd.to_datetime(establishments_df['opening_date']) <= pd.to_datetime(treatment_end))
        ]
        logger.info(f"✓ Filtered treatments before {treatment_end}")
    
    # Separate treatment sites from outcomes
    treatment_sites = establishments_df[establishments_df['type'] == treatment_type].copy()
    
    # Filter outcome types
    if outcome_types:
        outcomes = establishments_df[establishments_df['type'].isin(outcome_types)].copy()
    else:
        outcomes = establishments_df[establishments_df['type'] != treatment_type].copy()
    
    # Apply category filter
    if category_filter:
        filtered_outcomes = []
        for est_type, categories in category_filter.items():
            filtered = outcomes[
                (outcomes['type'] == est_type) & 
                (outcomes['category'].isin(categories))
            ]
            filtered_outcomes.append(filtered)
        outcomes = pd.concat(filtered_outcomes, ignore_index=True)
        logger.info(f"✓ Filtered to specific categories")
    
    logger.info(f"  Treatment sites ({treatment_type}): {len(treatment_sites)}")
    for outcome_type in outcomes['type'].unique():
        count = len(outcomes[outcomes['type'] == outcome_type])
        logger.info(f"  Outcomes ({outcome_type}): {count}")
    
    if len(treatment_sites) == 0:
        print(f"\n{Colors.RED}✗ No treatment sites found in data{Colors.ENDC}")
        return
    
    # Create GeoDataFrames
    print(f"\n{Colors.BOLD}[2/5] Creating spatial geometries...{Colors.ENDC}")
    treatment_gdf = gpd.GeoDataFrame(
        treatment_sites,
        geometry=[Point(lon, lat) for lon, lat in zip(treatment_sites['longitude'], treatment_sites['latitude'])],
        crs='EPSG:4326'
    )
    
    outcomes_gdf = gpd.GeoDataFrame(
        outcomes,
        geometry=[Point(lon, lat) for lon, lat in zip(outcomes['longitude'], outcomes['latitude'])],
        crs='EPSG:4326'
    )
    logger.info("✓ Created geometries")
    
    # Create one AOT around each treatment site
    print(f"\n{Colors.BOLD}[3/5] Creating Areas of Treatment (AOTs)...{Colors.ENDC}")
    aots_gdf = create_aots(treatment_gdf, aot_radius_m=aot_radius_m)
    
    # Spatial join
    print(f"\n{Colors.BOLD}[4/5] Performing spatial join...{Colors.ENDC}")
    spatial_join_df = spatial_join_establishments(aots_gdf, outcomes_gdf)
    
    # Filter sites with minimum establishments
    if min_establishments > 0:
        site_counts = spatial_join_df.groupby('site_id').size()
        valid_sites = site_counts[site_counts >= min_establishments].index
        spatial_join_df = spatial_join_df[spatial_join_df['site_id'].isin(valid_sites)]
        logger.info(f"✓ Filtered to sites with ≥{min_establishments} establishments: {len(valid_sites)} sites")
    
    # Generate panel
    print(f"\n{Colors.BOLD}[5/5] Generating panel dataset...{Colors.ENDC}")
    panel_df = generate_panel(spatial_join_df, start_date, end_date)
    
    # Create output directory
    timestamp = datetime.now().strftime('%d%m%y_%H%M%S')
    output_dir_path = Path(output_dir) / 'panels' / f'panel_{timestamp}'
    output_dir_path.mkdir(parents=True, exist_ok=True)
    
    # Save panel CSV
    panel_file = output_dir_path / 'panel.csv'
    panel_df.to_csv(panel_file, index=False)
    logger.info(f"✓ Saved panel to {panel_file}")
    
    # Save metadata
    metadata_file = output_dir_path / 'metadata.txt'
    with open(metadata_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("OSM PANEL DATASET METADATA\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Generation Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Config File: {config_path}\n")
        f.write(f"Input File: {input_file}\n\n")
        
        f.write("PANEL STRUCTURE:\n")
        f.write(f"  Total Observations: {len(panel_df):,}\n")
        f.write(f"  Unique Sites: {panel_df['site_id'].nunique()}\n")
        f.write(f"  AOT Radius: {aot_radius_m}m\n")
        f.write(f"  Time Periods: {panel_df['month'].nunique()} months\n")
        f.write(f"  Date Range: {start_date} to {end_date}\n\n")
        
        f.write("TREATMENT:\n")
        f.write(f"  Treatment Type: {treatment_type} openings\n")
        f.write(f"  Number of Sites: {len(treatment_sites)}\n")
        f.write(f"  Number of Cohorts: {panel_df['cohort_id'].nunique()}\n")
        f.write(f"  Cohorts: {', '.join(sorted(panel_df['cohort_id'].unique()))}\n")
        if treatment_start:
            f.write(f"  Treatment Start Filter: {treatment_start}\n")
        if treatment_end:
            f.write(f"  Treatment End Filter: {treatment_end}\n")
        f.write("\n")
        
        f.write("OUTCOMES:\n")
        outcome_types_str = ', '.join(sorted(outcomes['type'].unique()))
        f.write(f"  Establishment Types: {outcome_types_str}\n")
        f.write(f"  Measurements: Count, Entries, Exits\n")
        if category_filter:
            f.write(f"  Category Filter Applied: Yes\n")
            for est_type, cats in category_filter.items():
                f.write(f"    {est_type}: {', '.join(cats)}\n")
        if min_establishments > 0:
            f.write(f"  Minimum Establishments per Site: {min_establishments}\n")
        f.write("\n")
        
        f.write("PANEL VARIABLES:\n")
        f.write("  - site_id: Alphanumeric identifier for each treatment site\n")
        f.write("  - month: Month and year of observation (YYYY-MM)\n")
        f.write("  - poi_poc: Type of establishment (POI or POC)\n")
        f.write("  - count: Number of active establishments in month\n")
        f.write("  - entries: Number of new establishments in month\n")
        f.write("  - exits: Number of closed establishments in month\n")
        f.write("  - date_of_treatment: Month and year of treatment\n")
        f.write("  - cohort_id: Year and quarter of treatment (YYYYQX)\n")
        f.write("  - time_to_treatment: Months to/from treatment (k=0 is treatment month)\n\n")
        
        f.write("ANALYSIS READY:\n")
        f.write("  ✓ Long POI/POC panel\n")
        f.write("  ✓ Staggered treatment timing\n")
        f.write("  ✓ Event-study estimation\n")
        f.write("  ✓ PPML-compatible outcomes\n\n")
        
        f.write("="*80 + "\n")
    
    logger.info(f"✓ Saved metadata to {metadata_file}")
    
    # Save config copy
    if config:
        config_copy = output_dir_path / 'config.yml'
        with open(config_copy, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        logger.info(f"✓ Saved config copy to {config_copy}")
    
    # Print summary
    print("\n" + "="*80)
    print(f"{Colors.BOLD}{Colors.GREEN}✓ PANEL GENERATION COMPLETE{Colors.ENDC}")
    print("="*80)
    print(f"\n{Colors.BOLD}Panel Summary:{Colors.ENDC}")
    print(f"  Observations: {Colors.CYAN}{len(panel_df):,}{Colors.ENDC}")
    print(f"  Sites: {Colors.CYAN}{panel_df['site_id'].nunique()}{Colors.ENDC}")
    print(f"  Date range: {Colors.CYAN}{panel_df['month'].min()}{Colors.ENDC} to {Colors.CYAN}{panel_df['month'].max()}{Colors.ENDC}")
    print(f"  Cohorts: {Colors.CYAN}{panel_df['cohort_id'].nunique()}{Colors.ENDC}")
    
    print(f"\n{Colors.BOLD}Output:{Colors.ENDC} {Colors.CYAN}{output_dir_path}{Colors.ENDC}")
    print(f"  {Colors.GREEN}→{Colors.ENDC} panel.csv")
    print(f"  {Colors.GREEN}→{Colors.ENDC} metadata.txt")
    print(f"  {Colors.GREEN}→{Colors.ENDC} config.yml")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
