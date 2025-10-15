#!/usr/bin/env python3
# -------------------------------------------------------------------------------------------------
# DDD Panel Data Processor
#
# This script processes raw DDD sample data (CSV + GPKG) and creates a structured panel dataset
# for Difference-in-Difference-in-Differences analysis.
#
# Directory structure:
#   data/raw/ddd_sample_*/        - Raw outputs (CSV + GPKG)
#   data/processed/ddd_panel_*.csv - Processed panel data
#
# CLI usage:
#   python python/process_ddd_panel.py --input data/raw/ddd_sample_TIMESTAMP_RADIUS
#
# Requirements:
#   pip install pandas numpy geopandas shapely scikit-learn
# -------------------------------------------------------------------------------------------------

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from sklearn.cluster import DBSCAN

warnings.filterwarnings('ignore')

# ================================ CONFIGURATION ================================

# Panel time bounds (inclusive)
PANEL_START = "2016-01-01"
PANEL_END = "2025-09-01"

# Category mappings (OSM tags → canonical categories)
CATEGORY_MAPPING = {
    'bakery': 'bakery',
    'butcher': 'butcher',
    'greengrocer': 'greengrocer',
    'convenience': 'convenience',
    'funeral_directors': 'funeral_home',
    'hairdresser': 'hairdresser',
    'restaurant': 'restaurant',
}

# Treatment assignment
AFFECTED_CATEGORIES = {'bakery', 'butcher', 'greengrocer', 'convenience'}
PLACEBO_CATEGORIES = {'funeral_home', 'hairdresser', 'restaurant'}

# Ring definitions (closed-on-right convention)
RINGS = {
    'in': (0, 500),      # [0, 500] m
    'mid': (500, 1000),  # (500, 1000] m
    'out': (1000, 1500), # (1000, 1500] m
}

# Site clustering distance (meters)
SITE_CLUSTER_DISTANCE = 75

# Eligibility: minimum months pre and post
MIN_PRE_MONTHS = 12
MIN_POST_MONTHS = 12

# Event time caps
EVENT_TIME_MIN = -24
EVENT_TIME_MAX = 24

# ================================ UTILITIES ================================

def info(msg: str):
    print(f"[INFO] {msg}")

def warn(msg: str):
    print(f"[WARN] {msg}")

def error(msg: str):
    print(f"[ERROR] {msg}")

def success(msg: str):
    print(f"[OK] {msg}")

# ================================ SITE CLUSTERING ================================

def cluster_supermarkets(df_super: pd.DataFrame, eps_m: float = SITE_CLUSTER_DISTANCE) -> pd.DataFrame:
    """
    Cluster supermarkets within eps_m meters and assign site_id to each cluster.
    Returns df_super with added 'site_id' column.
    """
    info(f"Clustering {len(df_super)} supermarkets with eps={eps_m}m...")
    
    if len(df_super) == 0:
        return df_super.assign(site_id=pd.Series(dtype=str))
    
    # Use UTM coordinates for metric clustering
    coords = df_super[['lon', 'lat']].values
    
    # Simple geographic distance approximation (for small areas)
    # Convert to approximate meters: 1 degree ≈ 111km at equator
    lat_mean = coords[:, 1].mean()
    lon_scale = 111000 * np.cos(np.radians(lat_mean))
    lat_scale = 111000
    
    coords_m = np.column_stack([
        coords[:, 0] * lon_scale,
        coords[:, 1] * lat_scale
    ])
    
    # DBSCAN clustering
    clustering = DBSCAN(eps=eps_m, min_samples=1, metric='euclidean')
    labels = clustering.fit_predict(coords_m)
    
    # Create site IDs
    site_ids = [f"S{label:04d}" for label in labels]
    df_super = df_super.copy()
    df_super['site_id'] = site_ids
    
    n_sites = len(set(site_ids))
    info(f"Created {n_sites} supermarket sites from {len(df_super)} stores")
    
    return df_super

# ================================ SITE ENTRY DATES ================================

def compute_site_entry_dates(df_super: pd.DataFrame) -> pd.DataFrame:
    """
    For each site_id, compute t0 = earliest timestamp across all supermarkets in the cluster.
    Returns DataFrame with columns: site_id, t0
    """
    info("Computing site entry dates (t0)...")
    
    df_super['timestamp_dt'] = pd.to_datetime(df_super['timestamp'])
    
    site_entry = df_super.groupby('site_id')['timestamp_dt'].min().reset_index()
    site_entry.columns = ['site_id', 't0']
    
    # Convert to first day of month
    site_entry['t0'] = site_entry['t0'].dt.to_period('M').dt.to_timestamp()
    
    info(f"Computed entry dates for {len(site_entry)} sites")
    info(f"Entry date range: {site_entry['t0'].min()} to {site_entry['t0'].max()}")
    
    return site_entry

# ================================ NEAREST SITE ASSIGNMENT ================================

def assign_to_nearest_site(df_pois: pd.DataFrame, df_sites: pd.DataFrame, 
                          max_distance_m: float = 1500) -> pd.DataFrame:
    """
    Assign each POI/POC to the nearest supermarket site within max_distance_m.
    Returns df_pois with added columns: site_id, distance_to_site_m, ring
    """
    info(f"Assigning {len(df_pois)} POIs/POCs to nearest sites (max {max_distance_m}m)...")
    
    results = []
    
    # Get representative point for each site (centroid of cluster)
    site_coords = df_sites.groupby('site_id')[['lon', 'lat']].mean().reset_index()
    
    lat_mean = df_pois['lat'].mean()
    lon_scale = 111000 * np.cos(np.radians(lat_mean))
    lat_scale = 111000
    
    for idx, poi in df_pois.iterrows():
        poi_lon, poi_lat = poi['lon'], poi['lat']
        
        # Calculate distances to all sites
        dx = (site_coords['lon'] - poi_lon) * lon_scale
        dy = (site_coords['lat'] - poi_lat) * lat_scale
        distances = np.sqrt(dx**2 + dy**2)
        
        min_dist_idx = distances.idxmin()
        min_dist = distances[min_dist_idx]
        
        if min_dist <= max_distance_m:
            nearest_site = site_coords.loc[min_dist_idx, 'site_id']
            
            # Determine ring (closed-on-right)
            if min_dist <= 500:
                ring = 'in'
            elif min_dist <= 1000:
                ring = 'mid'
            else:
                ring = 'out'
            
            results.append({
                'osm_id': poi['osm_id'],
                'site_id': nearest_site,
                'distance_to_site_m': min_dist,
                'ring': ring
            })
    
    df_assigned = pd.DataFrame(results)
    
    if len(df_assigned) == 0:
        warn("No POIs/POCs assigned to any site!")
        return df_pois.assign(site_id=pd.Series(dtype=str),
                             distance_to_site_m=pd.Series(dtype=float),
                             ring=pd.Series(dtype=str))
    
    # Merge back with original POI data
    df_result = df_pois.merge(df_assigned, on='osm_id', how='inner')
    
    info(f"Assigned {len(df_result)} POIs/POCs to sites (dropped {len(df_pois) - len(df_result)} beyond {max_distance_m}m)")
    
    # Ring distribution
    ring_counts = df_result['ring'].value_counts().to_dict()
    for ring, count in ring_counts.items():
        info(f"  Ring '{ring}': {count} POIs/POCs")
    
    return df_result

# ================================ PANEL CONSTRUCTION ================================

def create_area_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Add area_id = site_id + '_' + ring"""
    df['area_id'] = df['site_id'] + '_' + df['ring']
    return df

def build_balanced_panel(df_assigned: pd.DataFrame, df_site_entry: pd.DataFrame,
                        panel_start: str, panel_end: str) -> pd.DataFrame:
    """
    Build balanced panel with all area_id × category × month combinations.
    """
    info("Building balanced panel...")
    
    # Standardize categories
    df_assigned['category'] = df_assigned['subcategory'].map(CATEGORY_MAPPING)
    df_assigned = df_assigned[df_assigned['category'].notna()]
    
    # Create area_ids
    df_assigned = create_area_ids(df_assigned)
    
    # Get first appearance month for each POI
    df_assigned['first_month'] = pd.to_datetime(df_assigned['timestamp']).dt.to_period('M').dt.to_timestamp()
    
    # Time grid
    time_grid = pd.date_range(start=panel_start, end=panel_end, freq='MS')
    
    # Get all area_ids
    area_info = df_assigned[['area_id', 'site_id', 'ring']].drop_duplicates()
    
    # Get all categories
    all_categories = list(set(CATEGORY_MAPPING.values()))
    
    info(f"Panel dimensions: {len(area_info)} areas × {len(all_categories)} categories × {len(time_grid)} months")
    
    # Create full panel skeleton
    panel_rows = []
    
    for _, area in area_info.iterrows():
        area_id = area['area_id']
        site_id = area['site_id']
        ring = area['ring']
        
        for category in all_categories:
            for month in time_grid:
                panel_rows.append({
                    'area_id': area_id,
                    'site_id': site_id,
                    'ring': ring,
                    'category': category,
                    'month': month
                })
    
    panel_df = pd.DataFrame(panel_rows)
    info(f"Created panel skeleton: {len(panel_df):,} rows")
    
    # Compute stock counts (y)
    info("Computing stock counts (y)...")
    
    # For each area_id × category, get list of POIs with their first_month
    poi_stocks = df_assigned.groupby(['area_id', 'category', 'osm_id'])['first_month'].min().reset_index()
    
    # For each panel row, count POIs that have appeared by that month
    stock_counts = []
    
    for area_id in panel_df['area_id'].unique():
        for category in panel_df['category'].unique():
            # Get POIs in this area × category
            mask = (poi_stocks['area_id'] == area_id) & (poi_stocks['category'] == category)
            pois_here = poi_stocks[mask]
            
            if len(pois_here) == 0:
                continue
            
            # For each month, count how many POIs have appeared
            for month in time_grid:
                count = (pois_here['first_month'] <= month).sum()
                if count > 0:
                    stock_counts.append({
                        'area_id': area_id,
                        'category': category,
                        'month': month,
                        'y': count
                    })
    
    stock_df = pd.DataFrame(stock_counts)
    
    # Merge with panel skeleton
    panel_df = panel_df.merge(stock_df, on=['area_id', 'category', 'month'], how='left')
    panel_df['y'] = panel_df['y'].fillna(0).astype(int)
    
    info(f"Non-zero observations: {(panel_df['y'] > 0).sum():,} / {len(panel_df):,}")
    
    # Add site entry dates (t0)
    panel_df = panel_df.merge(df_site_entry, on='site_id', how='left')
    
    # Add treatment indicator (A)
    panel_df['A'] = panel_df['category'].apply(lambda c: 1 if c in AFFECTED_CATEGORIES else 0)
    
    # Add ring metadata
    panel_df['ring_min_m'] = panel_df['ring'].map({'in': 0, 'mid': 500, 'out': 1000})
    panel_df['ring_max_m'] = panel_df['ring'].map({'in': 500, 'mid': 1000, 'out': 1500})
    
    # Compute exposure (E)
    panel_df['E'] = (panel_df['month'] >= panel_df['t0']).astype(int)
    
    # Compute event_time
    panel_df['event_time'] = (
        12 * (panel_df['month'].dt.year - panel_df['t0'].dt.year) +
        (panel_df['month'].dt.month - panel_df['t0'].dt.month)
    )
    
    # Compute ev_cap
    panel_df['ev_cap'] = panel_df['event_time'].clip(EVENT_TIME_MIN, EVENT_TIME_MAX)
    
    # Add panel metadata
    panel_df['sample_window_start'] = pd.to_datetime(panel_start)
    panel_df['sample_window_end'] = pd.to_datetime(panel_end)
    
    # Compute eligibility
    panel_start_dt = pd.to_datetime(panel_start)
    panel_end_dt = pd.to_datetime(panel_end)
    
    def is_eligible(t0):
        if pd.isna(t0):
            return 0
        pre_months = (t0 - panel_start_dt).days / 30.44
        post_months = (panel_end_dt - t0).days / 30.44
        return 1 if (pre_months >= MIN_PRE_MONTHS and post_months >= MIN_POST_MONTHS) else 0
    
    # Eligibility is site-level
    site_eligibility = panel_df.groupby('site_id')['t0'].first().apply(is_eligible)
    site_eligibility_df = site_eligibility.reset_index()
    site_eligibility_df.columns = ['site_id', 'eligible']
    
    panel_df = panel_df.merge(site_eligibility_df, on='site_id', how='left')
    
    eligible_count = panel_df['eligible'].sum() / len(panel_df['month'].unique()) / len(all_categories)
    info(f"Eligible sites: {int(eligible_count)} / {len(area_info['site_id'].unique())}")
    
    return panel_df

# ================================ AREA CALCULATION ================================

def compute_ring_areas(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute approximate geodesic area for each ring.
    For simplicity, using circular approximation.
    """
    info("Computing ring areas...")
    
    # Circular ring areas (km²)
    ring_areas = {
        'in': np.pi * (0.5**2),              # Circle of 500m radius
        'mid': np.pi * (1.0**2 - 0.5**2),    # Annulus 500-1000m
        'out': np.pi * (1.5**2 - 1.0**2),    # Annulus 1000-1500m
    }
    
    df['area_km2'] = df['ring'].map(ring_areas)
    
    return df

# ================================ MAIN PROCESSING ================================

def process_ddd_sample(input_dir: str, output_dir: str):
    """Main processing pipeline"""
    
    print("="*70)
    print("DDD PANEL DATA PROCESSOR")
    print("="*70)
    
    input_path = Path(input_dir)
    
    if not input_path.exists():
        error(f"Input directory not found: {input_dir}")
        sys.exit(1)
    
    # Read raw data
    csv_file = input_path / "ddd_sample.csv"
    gpkg_file = input_path / "ddd_sample.gpkg"
    
    if not csv_file.exists():
        error(f"CSV file not found: {csv_file}")
        sys.exit(1)
    
    info(f"Reading raw data from: {input_path}")
    
    # Read CSV (skip comment lines)
    df_raw = pd.read_csv(csv_file, comment='#')
    info(f"Loaded {len(df_raw)} raw observations")
    
    # Split by kind
    df_supermarkets = df_raw[df_raw['kind'] == 'supermarket'].copy()
    df_pois = df_raw[df_raw['kind'].isin(['poi', 'control'])].copy()
    
    info(f"Supermarkets: {len(df_supermarkets)}, POIs/POCs: {len(df_pois)}")
    
    # Step 1: Cluster supermarkets into sites
    df_supermarkets = cluster_supermarkets(df_supermarkets)
    
    # Step 2: Compute site entry dates
    df_site_entry = compute_site_entry_dates(df_supermarkets)
    
    # Step 3: Assign POIs to nearest site
    df_pois_assigned = assign_to_nearest_site(df_pois, df_supermarkets)
    
    if len(df_pois_assigned) == 0:
        error("No POIs assigned to sites. Cannot build panel.")
        sys.exit(1)
    
    # Step 4: Build balanced panel
    panel_df = build_balanced_panel(df_pois_assigned, df_site_entry, PANEL_START, PANEL_END)
    
    # Step 5: Compute ring areas
    panel_df = compute_ring_areas(panel_df)
    
    # Step 6: Validate panel
    info("Validating panel integrity...")
    
    # Check uniqueness
    duplicates = panel_df.duplicated(subset=['area_id', 'category', 'month']).sum()
    if duplicates > 0:
        warn(f"Found {duplicates} duplicate rows!")
    else:
        success("Panel uniqueness: OK")
    
    # Check completeness
    expected_months = len(pd.date_range(PANEL_START, PANEL_END, freq='MS'))
    actual_months = panel_df.groupby(['area_id', 'category']).size()
    incomplete = (actual_months != expected_months).sum()
    if incomplete > 0:
        warn(f"Found {incomplete} incomplete area×category series")
    else:
        success("Panel completeness: OK")
    
    # Check for NAs in critical columns
    critical_cols = ['area_id', 'site_id', 'month', 'category', 'A', 'ring', 'y', 'E', 't0', 'event_time', 'ev_cap']
    na_counts = panel_df[critical_cols].isna().sum()
    if na_counts.sum() > 0:
        warn(f"Missing values found:\n{na_counts[na_counts > 0]}")
    else:
        success("No missing values in critical columns")
    
    # Step 7: Reorder and select final columns
    final_columns = [
        # Primary keys
        'area_id', 'site_id', 'month', 'category',
        # Treatment flags
        'A', 't0', 'E', 'event_time', 'ev_cap',
        # Spatial design
        'ring', 'ring_min_m', 'ring_max_m', 'area_km2',
        # Outcome
        'y',
        # Eligibility
        'eligible',
        # Metadata
        'sample_window_start', 'sample_window_end'
    ]
    
    panel_final = panel_df[final_columns].copy()
    
    # Sort by area_id, category, month
    panel_final = panel_final.sort_values(['area_id', 'category', 'month']).reset_index(drop=True)
    
    # Step 8: Save processed data
    os.makedirs(output_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    radius = input_path.name.split('_')[-1]
    output_file = os.path.join(output_dir, f"ddd_panel_{timestamp}_{radius}.csv")
    
    panel_final.to_csv(output_file, index=False)
    success(f"Panel data saved: {output_file}")
    
    # Summary statistics
    print("\n" + "="*70)
    print("PANEL SUMMARY")
    print("="*70)
    info(f"Total rows: {len(panel_final):,}")
    info(f"Unique areas: {panel_final['area_id'].nunique()}")
    info(f"Unique sites: {panel_final['site_id'].nunique()}")
    info(f"Categories: {panel_final['category'].nunique()}")
    info(f"Time periods: {panel_final['month'].nunique()} months")
    info(f"Date range: {panel_final['month'].min()} to {panel_final['month'].max()}")
    info(f"Affected observations (A=1): {(panel_final['A'] == 1).sum():,}")
    info(f"Placebo observations (A=0): {(panel_final['A'] == 0).sum():,}")
    info(f"Exposed observations (E=1): {(panel_final['E'] == 1).sum():,}")
    info(f"Eligible observations: {(panel_final['eligible'] == 1).sum():,}")
    info(f"Non-zero outcomes: {(panel_final['y'] > 0).sum():,} ({100*(panel_final['y'] > 0).mean():.1f}%)")
    
    print("\nRing distribution:")
    for ring in ['in', 'mid', 'out']:
        count = (panel_final['ring'] == ring).sum() // panel_final['month'].nunique() // panel_final['category'].nunique()
        print(f"  {ring}: {count} areas")
    
    print("\nCategory distribution (affected A=1):")
    for cat in AFFECTED_CATEGORIES:
        mask = (panel_final['category'] == cat) & (panel_final['A'] == 1)
        print(f"  {cat}: {mask.sum():,} obs")
    
    print("\nCategory distribution (placebo A=0):")
    for cat in PLACEBO_CATEGORIES:
        mask = (panel_final['category'] == cat) & (panel_final['A'] == 0)
        print(f"  {cat}: {mask.sum():,} obs")
    
    print("="*70)
    success("Processing complete!")
    
    return output_file

# ================================ CLI ================================

def main():
    parser = argparse.ArgumentParser(description="Process DDD sample data into panel format")
    parser.add_argument("--input", required=True, help="Input directory (raw/ddd_sample_*)")
    parser.add_argument("--output", default="data/processed", help="Output directory for panel data")
    
    args = parser.parse_args()
    
    try:
        process_ddd_sample(args.input, args.output)
    except Exception as e:
        error(f"Processing failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
