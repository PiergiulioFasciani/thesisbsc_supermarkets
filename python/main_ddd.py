# -------------------------------------------------------------------------------------------------
# DDD (Difference-in-Difference-in-Differences) OSM Data Extraction
#
# This script extracts supermarkets and POIs from OSM data within an AOI and a buffer zone,
# along with their creation timestamps (as proxy for opening dates) for DDD analysis.
#
# Output: CSV file with all supermarkets and POIs, their locations, types, and timestamps
# Format: ddd_sample_DDMMYYYY_HHMMSS_[radiusAOI].csv
#
# CLI usage:
#   python python/main_did.py --center "45.43318,9.18378" --radius-m 1200 \
#       --pbf data/pbf/nord-ovest-latest.osm.pbf --out-dir data
#
# Requirements:
#   pip install pyosmium geopandas shapely pyproj pandas numpy tqdm rich
# -------------------------------------------------------------------------------------------------

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon
from shapely.ops import transform as shp_transform
from shapely import wkb as shp_wkb
from pyproj import CRS, Transformer

# Pretty console (optional)
try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
    RICH = True
    console = Console()
except Exception:
    RICH = False
    console = None

try:
    from tqdm import tqdm
    TQDM = True
except Exception:
    TQDM = False
    def tqdm(x, **kwargs): return x

# PyOsmium
try:
    import osmium as osm
except ImportError:
    raise SystemExit("pyosmium is required. Install with: pip install pyosmium")

# ------------------------------------ CONFIG / TAGS ------------------------------------

DUOMO_LAT, DUOMO_LON = 45.464211, 9.191383

# Target POIs (limited to specified categories for DDD analysis)
POI_SHOPS = {"bakery", "greengrocer", "butcher", "convenience"}

# Control POIs (for DDD analysis - only in AOI)
CONTROL_SHOPS = {"funeral_directors", "hairdresser"}  # funeral_directors = funeral homes/services
CONTROL_AMENITIES = {"restaurant"}

# Supermarkets (separate category)
SUPERMARKET_SHOPS = {"supermarket"}

# ---- Portable defaults ----
DEFAULT_PBF = os.path.join("data", "pbf", "nord-ovest-latest.osm.pbf")
DEFAULT_OUT_BASE = os.path.join("data", "raw")

# Buffer distance (1000m beyond AOI as specified)
BUFFER_DISTANCE = 1000

# ----------------------------------------- UX -----------------------------------------

def banner(title: str):
    if RICH: 
        console.print(Panel.fit(f"[bold cyan]{title}[/]", border_style="cyan"))
    else: 
        print(f"\n=== {title} ===")

def info(msg: str):
    if RICH: 
        console.print(f"[bold]• {msg}[/]")
    else: 
        print(f"• {msg}")

def success(msg: str):
    if RICH: 
        console.print(f"[green]✓ {msg}[/]")
    else: 
        print(f"[OK] {msg}")

def warn(msg: str):
    if RICH: 
        console.print(f"[yellow]! {msg}[/]")
    else: 
        print(f"[WARN] {msg}")

def err(msg: str):
    if RICH: 
        console.print(f"[red]✗ {msg}[/]")
    else: 
        print(f"[ERR] {msg}")

# ------------------------------------- HELPERS -------------------------------------

def parse_center(s: str) -> Tuple[float, float]:
    s = s.strip().replace(" ", "")
    lat_s, lon_s = s.split(",", 1)
    return float(lat_s), float(lon_s)

def make_aeqd_crs(center_lon: float, center_lat: float) -> CRS:
    proj4 = f"+proj=aeqd +lat_0={center_lat} +lon_0={center_lon} +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    return CRS.from_proj4(proj4)

def circle_polygon_wgs84(center_lon: float, center_lat: float, radius_m: float) -> Polygon:
    wgs84 = CRS.from_epsg(4326)
    aeqd = make_aeqd_crs(center_lon, center_lat)
    to_aeqd = Transformer.from_crs(wgs84, aeqd, always_xy=True).transform
    to_wgs = Transformer.from_crs(aeqd, wgs84, always_xy=True).transform
    poly_m = shp_transform(to_aeqd, Point(center_lon, center_lat)).buffer(radius_m, resolution=64)
    poly = shp_transform(to_wgs, poly_m)
    try: 
        poly = poly.buffer(0)
    except Exception: 
        pass
    return poly

def prepare_output_folder(base_dir: str, radius_m: float) -> Tuple[str, str, str]:
    """
    Create output folder and filenames like:
      ddd_sample_DDMMYYYY_HHMMSS_<AOI meters>/
        ddd_sample.csv
        ddd_sample.gpkg
    Example: ddd_sample_06102025_143918_1200/
    """
    ts = datetime.now().strftime("%d%m%Y_%H%M%S")
    radius_tag = f"{int(round(radius_m))}"
    folder_name = f"ddd_sample_{ts}_{radius_tag}"
    folder_path = os.path.join(base_dir, folder_name)
    
    os.makedirs(folder_path, exist_ok=True)
    
    csv_file = os.path.join(folder_path, "ddd_sample.csv")
    gpkg_file = os.path.join(folder_path, "ddd_sample.gpkg")
    
    return folder_path, csv_file, gpkg_file

def local_utm_epsg(lon: float, lat: float) -> int:
    """Get the appropriate UTM EPSG code for the given coordinates"""
    zone = int((lon + 180.0) / 6.0) + 1
    return (32600 if lat >= 0 else 32700) + zone

def calculate_distance_to_duomo(lon: float, lat: float, duomo_lon: float = DUOMO_LON, duomo_lat: float = DUOMO_LAT) -> float:
    """Calculate distance from point to Duomo in meters using UTM projection"""
    # Use UTM projection for accurate distance calculation
    utm_epsg = local_utm_epsg(lon, lat)
    crs_m = CRS.from_epsg(utm_epsg)
    
    to_m = Transformer.from_crs(4326, crs_m, always_xy=True).transform
    
    # Transform both points to UTM
    point_m = shp_transform(to_m, Point(lon, lat))
    duomo_m = shp_transform(to_m, Point(duomo_lon, duomo_lat))
    
    return float(point_m.distance(duomo_m))

# ------------------------------------- OSM PARSER -------------------------------------

@dataclass
class DDDFeature:
    id: int
    lon: float
    lat: float
    name: Optional[str]
    kind: str          # 'supermarket', 'poi'
    subcategory: str   # specific tag value (e.g., 'bakery', 'supermarket')
    source: str        # 'node' | 'area'
    timestamp: Optional[str]  # OSM timestamp as string
    zone: str         # 'aoi' or 'buffer'
    distance_to_duomo_m: float  # Distance to Duomo in meters

class DDDTagMatcher:
    def match_kind_and_subcategory(self, tags) -> Optional[Tuple[str, str]]:
        shop = tags.get("shop")
        amenity = tags.get("amenity")

        if shop in SUPERMARKET_SHOPS:
            return "supermarket", shop
        if shop in POI_SHOPS:
            return "poi", shop
        if shop in CONTROL_SHOPS:
            return "control", shop
        if amenity in CONTROL_AMENITIES:
            return "control", amenity
        
        return None

class DDDAreaFilter:
    def __init__(self, aoi_poly: Polygon, buffer_poly: Polygon):
        self.aoi = aoi_poly
        self.buffer = buffer_poly
        self.aoi_bounds = aoi_poly.bounds
        self.buffer_bounds = buffer_poly.bounds
    
    def classify_lonlat(self, lon: float, lat: float) -> Optional[str]:
        """Returns 'aoi', 'buffer', or None"""
        # Quick bounds check first
        minx, miny, maxx, maxy = self.buffer_bounds
        if not (minx <= lon <= maxx and miny <= lat <= maxy):
            return None
        
        point = Point(lon, lat)
        
        if self.aoi.covers(point):
            return "aoi"
        elif self.buffer.covers(point):
            return "buffer"
        else:
            return None

class DDDOSMCollector(osm.SimpleHandler):
    """ Collect supermarkets and POIs in AOI and buffer zone with timestamps """
    
    def __init__(self, tm: DDDTagMatcher, af: DDDAreaFilter):
        super().__init__()
        self.tm = tm
        self.af = af
        self.wkbf = osm.geom.WKBFactory()
        self.features: List[DDDFeature] = []
        
        # Cache buffer bounds for efficiency
        self.buffer_bounds = af.buffer_bounds

    def _format_timestamp(self, osm_obj) -> Optional[str]:
        """Extract and format OSM timestamp"""
        try:
            if hasattr(osm_obj, 'timestamp') and osm_obj.timestamp:
                return osm_obj.timestamp.strftime("%Y-%m-%d %H:%M:%S")
        except:
            pass
        return None

    def _maybe_keep(self, feature: DDDFeature):
        if feature.zone:  # zone is set if point is in AOI or buffer
            # Control features only in AOI (not buffer)
            if feature.kind == "control" and feature.zone != "aoi":
                return
            self.features.append(feature)

    def node(self, n):
        match_result = self.tm.match_kind_and_subcategory(n.tags)
        if match_result is None or not n.location.valid(): 
            return
        
        lon, lat = n.location.lon, n.location.lat
        
        # Early filtering - skip if outside buffer bounds completely
        minx, miny, maxx, maxy = self.buffer_bounds
        if not (minx <= lon <= maxx and miny <= lat <= maxy):
            return
        
        kind, subcategory = match_result
        zone = self.af.classify_lonlat(lon, lat)
        
        if zone is None:
            return
        
        name = n.tags.get("name")
        timestamp = self._format_timestamp(n)
        distance_to_duomo = calculate_distance_to_duomo(lon, lat)
        
        feature = DDDFeature(
            id=n.id,
            lon=lon,
            lat=lat,
            name=name,
            kind=kind,
            subcategory=subcategory,
            source="node",
            timestamp=timestamp,
            zone=zone,
            distance_to_duomo_m=distance_to_duomo
        )
        self._maybe_keep(feature)

    def area(self, a):
        match_result = self.tm.match_kind_and_subcategory(a.tags)
        if match_result is None:
            return
        
        kind, subcategory = match_result
        
        try:
            wkb = None
            try: 
                wkb = self.wkbf.create_multipolygon(a)
            except Exception: 
                pass
            if wkb is None:
                try: 
                    wkb = self.wkbf.create_polygon(a)
                except Exception: 
                    return
            
            geom = shp_wkb.loads(wkb, hex=True)
            c = geom.centroid
            lon, lat = c.x, c.y
            
            # Early filtering - skip if outside buffer bounds completely
            minx, miny, maxx, maxy = self.buffer_bounds
            if not (minx <= lon <= maxx and miny <= lat <= maxy):
                return
            
            zone = self.af.classify_lonlat(lon, lat)
            
            if zone is None:
                return
            
            name = a.tags.get("name")
            timestamp = self._format_timestamp(a)
            distance_to_duomo = calculate_distance_to_duomo(lon, lat)
            
            feature = DDDFeature(
                id=a.id,
                lon=lon,
                lat=lat,
                name=name,
                kind=kind,
                subcategory=subcategory,
                source="area",
                timestamp=timestamp,
                zone=zone,
                distance_to_duomo_m=distance_to_duomo
            )
            self._maybe_keep(feature)
        except Exception:
            return

def parse_pbf_for_ddd(pbf_path: str, aoi_poly: Polygon, buffer_poly: Polygon):
    """Parse PBF file to extract supermarkets and POIs for DDD analysis"""
    banner("Step 2 — Reading OSM for DDD analysis (AOI + Buffer zone optimized)")
    tm = DDDTagMatcher()
    af = DDDAreaFilter(aoi_poly, buffer_poly)
    
    # Get bounding box of buffer zone for efficient pre-filtering
    minx, miny, maxx, maxy = buffer_poly.bounds
    info(f"Target area: lon [{minx:.6f}, {maxx:.6f}], lat [{miny:.6f}, {maxy:.6f}]")
    
    if RICH:
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task = progress.add_task("Parsing PBF for DDD…", total=None)
            h = DDDOSMCollector(tm, af)
            h.apply_file(str(pbf_path), locations=True, idx="flex_mem")
            progress.update(task, description="Finalizing…")
    else:
        h = DDDOSMCollector(tm, af)
        h.apply_file(str(pbf_path), locations=True, idx="flex_mem")

    info(f"Collected {len(h.features)} features for DDD analysis")
    
    # Count by zone and type
    aoi_count = sum(1 for f in h.features if f.zone == "aoi")
    buffer_count = sum(1 for f in h.features if f.zone == "buffer")
    supermarket_count = sum(1 for f in h.features if f.kind == "supermarket")
    poi_count = sum(1 for f in h.features if f.kind == "poi")
    control_count = sum(1 for f in h.features if f.kind == "control")
    
    info(f"AOI: {aoi_count}, Buffer: {buffer_count}")
    info(f"Supermarkets: {supermarket_count}, POIs: {poi_count}, Controls: {control_count}")
    
    return h.features

def features_to_dataframe(features: List[DDDFeature]) -> pd.DataFrame:
    """Convert features to pandas DataFrame"""
    if not features:
        return pd.DataFrame(columns=[
            'osm_id', 'name', 'kind', 'subcategory', 'source', 'zone',
            'lat', 'lon', 'timestamp', 'distance_to_duomo_m'
        ])
    
    data = []
    for f in features:
        data.append({
            'osm_id': f.id,
            'name': f.name,
            'kind': f.kind,
            'subcategory': f.subcategory,
            'source': f.source,
            'zone': f.zone,
            'lat': f.lat,
            'lon': f.lon,
            'timestamp': f.timestamp,
            'distance_to_duomo_m': f.distance_to_duomo_m
        })
    
    return pd.DataFrame(data)

def create_gpkg_output(features: List[DDDFeature], gpkg_path: str, 
                       aoi_poly: Polygon, buffer_poly: Polygon,
                       center_lon: float, center_lat: float):
    """Create GeoPackage with separate layers for different feature types"""
    banner("Step 4 — Creating GeoPackage output")
    
    # Use appropriate UTM projection for the area
    utm_epsg = local_utm_epsg(center_lon, center_lat)
    crs_m = CRS.from_epsg(utm_epsg)
    to_m = Transformer.from_crs(4326, crs_m, always_xy=True).transform
    
    # Convert features to GeoDataFrames by type
    pois_data = []
    pocs_data = []
    supermarkets_data = []
    
    for f in features:
        geom = Point(f.lon, f.lat)
        record = {
            'osm_id': f.id,
            'name': f.name,
            'subcategory': f.subcategory,
            'source': f.source,
            'zone': f.zone,
            'timestamp': f.timestamp,
            'distance_to_duomo_m': f.distance_to_duomo_m,
            'geometry': geom
        }
        
        if f.kind == 'poi':
            pois_data.append(record)
        elif f.kind == 'control':
            pocs_data.append(record)
        elif f.kind == 'supermarket':
            supermarkets_data.append(record)
    
    # Create GeoDataFrames
    gdf_pois = gpd.GeoDataFrame(pois_data, crs="EPSG:4326") if pois_data else gpd.GeoDataFrame(columns=['osm_id', 'name', 'subcategory', 'source', 'zone', 'timestamp', 'distance_to_duomo_m', 'geometry'], crs="EPSG:4326")
    gdf_pocs = gpd.GeoDataFrame(pocs_data, crs="EPSG:4326") if pocs_data else gpd.GeoDataFrame(columns=['osm_id', 'name', 'subcategory', 'source', 'zone', 'timestamp', 'distance_to_duomo_m', 'geometry'], crs="EPSG:4326")
    gdf_supermarkets = gpd.GeoDataFrame(supermarkets_data, crs="EPSG:4326") if supermarkets_data else gpd.GeoDataFrame(columns=['osm_id', 'name', 'subcategory', 'source', 'zone', 'timestamp', 'distance_to_duomo_m', 'geometry'], crs="EPSG:4326")
    
    # Project to UTM for consistent display
    gdf_pois = gdf_pois.to_crs(crs_m)
    gdf_pocs = gdf_pocs.to_crs(crs_m)
    gdf_supermarkets = gdf_supermarkets.to_crs(crs_m)
    
    # Create Duomo layer (marked as a star)
    duomo_point = Point(DUOMO_LON, DUOMO_LAT)
    gdf_duomo = gpd.GeoDataFrame(
        [{'name': 'Duomo di Milano', 'marker': 'star', 'geometry': duomo_point}],
        crs="EPSG:4326"
    ).to_crs(crs_m)
    
    # Create AOI and Buffer circle layers
    gdf_aoi = gpd.GeoDataFrame(
        [{'name': 'AOI', 'type': 'area_of_interest', 'geometry': aoi_poly}],
        crs="EPSG:4326"
    ).to_crs(crs_m)
    
    gdf_buffer = gpd.GeoDataFrame(
        [{'name': 'Buffer Zone', 'type': 'buffer', 'geometry': buffer_poly}],
        crs="EPSG:4326"
    ).to_crs(crs_m)
    
    # Write all layers to GeoPackage
    layer_info = []
    
    if len(gdf_pois) > 0:
        gdf_pois.to_file(gpkg_path, layer='pois', driver='GPKG')
        layer_info.append(f"POIs: {len(gdf_pois)}")
    
    if len(gdf_pocs) > 0:
        gdf_pocs.to_file(gpkg_path, layer='pocs', driver='GPKG', mode='a')
        layer_info.append(f"POCs: {len(gdf_pocs)}")
    
    if len(gdf_supermarkets) > 0:
        gdf_supermarkets.to_file(gpkg_path, layer='supermarkets', driver='GPKG', mode='a')
        layer_info.append(f"Supermarkets: {len(gdf_supermarkets)}")
    
    gdf_duomo.to_file(gpkg_path, layer='duomo', driver='GPKG', mode='a')
    layer_info.append("Duomo: 1 (star)")
    
    gdf_aoi.to_file(gpkg_path, layer='aoi_circle', driver='GPKG', mode='a')
    layer_info.append("AOI circle: 1")
    
    gdf_buffer.to_file(gpkg_path, layer='buffer_circle', driver='GPKG', mode='a')
    layer_info.append("Buffer circle: 1")
    
    success(f"GeoPackage written: {gpkg_path}")
    for info_line in layer_info:
        info(f"  • {info_line}")

# ------------------------------------- MAIN -------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Extract supermarkets and POIs for DDD analysis")
    ap.add_argument("--center", required=True, help="lat,lon (e.g. '45.43318,9.18378')")
    ap.add_argument("--radius-m", type=float, required=True, help="AOI radius in meters")
    ap.add_argument("--pbf", default=DEFAULT_PBF, help="Path to OSM PBF file")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_BASE, help="Output directory")
    args = ap.parse_args()

    # Step 1 — Setup
    banner("Step 1 — Setup for DDD Analysis")
    center_lat, center_lon = parse_center(args.center)
    pbf_path = Path(args.pbf).resolve()
    
    if not pbf_path.exists():
        err(f"PBF not found: {pbf_path}")
        raise SystemExit(1)

    # Create output directory and prepare filenames
    os.makedirs(args.out_dir, exist_ok=True)
    
    output_folder, csv_file, gpkg_file = prepare_output_folder(args.out_dir, args.radius_m)
    
    info(f"Center: ({center_lat:.6f}, {center_lon:.6f})  Radius: {args.radius_m:.0f} m")
    info(f"Buffer: {BUFFER_DISTANCE} m beyond AOI")
    info(f"PBF   : {pbf_path}")
    info(f"Output folder: {output_folder}")
    info(f"  • CSV:  {os.path.basename(csv_file)}")
    info(f"  • GPKG: {os.path.basename(gpkg_file)}")

    # Build AOI and buffer polygons first (WGS84)
    banner("Step 1a — Building AOI and Buffer zones")
    aoi_poly = circle_polygon_wgs84(center_lon=center_lon, center_lat=center_lat, radius_m=args.radius_m)
    buffer_poly = circle_polygon_wgs84(center_lon=center_lon, center_lat=center_lat, 
                                      radius_m=args.radius_m + BUFFER_DISTANCE)
    
    aoi_minx, aoi_miny, aoi_maxx, aoi_maxy = aoi_poly.bounds
    buf_minx, buf_miny, buf_maxx, buf_maxy = buffer_poly.bounds
    
    info(f"AOI bbox (WGS84): lon [{aoi_minx:.6f}, {aoi_maxx:.6f}], lat [{aoi_miny:.6f}, {aoi_maxy:.6f}]")
    info(f"Buffer bbox (WGS84): lon [{buf_minx:.6f}, {buf_maxx:.6f}], lat [{buf_miny:.6f}, {buf_maxy:.6f}]")
    
    # Calculate area reduction for efficiency info
    total_area_km2 = (buf_maxx - buf_minx) * (buf_maxy - buf_miny) * 111 * 111  # rough km²
    info(f"Parse area: ~{total_area_km2:.1f} km² (reduced from full PBF)")

    # Step 2 — Parse OSM
    features = parse_pbf_for_ddd(str(pbf_path), aoi_poly, buffer_poly)

    # Step 3 — Create DataFrame and export to CSV
    banner("Step 3 — Creating CSV output")
    df = features_to_dataframe(features)
    
    if len(df) == 0:
        warn("No features found! Creating empty CSV.")
    
    # Add some metadata as comments at the top of the CSV
    comments = [
        f"# DDD Analysis Dataset Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
        f"# AOI Center: {center_lat:.6f}, {center_lon:.6f}\n",
        f"# AOI Radius: {args.radius_m:.0f} meters\n",
        f"# Buffer Distance: {BUFFER_DISTANCE} meters\n",
        f"# Total Features: {len(df)}\n",
        f"# Source PBF: {pbf_path.name}\n",
        "#\n",
        "# Column descriptions:\n",
        "# - osm_id: OpenStreetMap element ID\n",
        "# - name: Name of the facility (if available)\n",
        "# - kind: Category (supermarket, poi, or control)\n",
        "# - subcategory: Specific OSM tag value (supermarket, bakery, butcher, etc.)\n",
        "# - source: OSM element type (node or area)\n",
        "# - zone: Whether in AOI or buffer zone\n",
        "# - lat, lon: Coordinates in WGS84\n",
        "# - timestamp: OSM creation timestamp (proxy for opening date)\n",
        "# - distance_to_duomo_m: Distance to Duomo di Milano in meters\n",
        "#\n",
        "# Treatment POI Categories: bakery, greengrocer, butcher, convenience\n",
        "# Supermarket Category: supermarket\n",
        "# Control Categories (AOI only): funeral_directors, hairdresser, restaurant\n",
        "#\n"
    ]
    
    with open(csv_file, "w", encoding="utf-8", newline="") as f:
        f.writelines(comments)
        df.to_csv(f, index=False, encoding="utf-8")
    
    success(f"CSV written: {csv_file}")
    
    # Step 4 — Create GeoPackage
    create_gpkg_output(features, gpkg_file, aoi_poly, buffer_poly, center_lon, center_lat)
    
    # Summary statistics
    if len(df) > 0:
        banner("Step 5 — Summary Statistics")
        info(f"Total features: {len(df)}")
        
        # By zone
        zone_counts = df['zone'].value_counts()
        for zone, count in zone_counts.items():
            info(f"  {zone.upper()}: {count}")
        
        # By kind
        kind_counts = df['kind'].value_counts()
        for kind, count in kind_counts.items():
            info(f"  {kind}: {count}")
        
        # Control subcategories breakdown
        if 'control' in kind_counts:
            control_subcats = df[df['kind'] == 'control']['subcategory'].value_counts()
            for subcat, count in control_subcats.items():
                info(f"    - {subcat}: {count}")
        
        # Distance to Duomo statistics
        poi_distances = df[df['kind'] == 'poi']['distance_to_duomo_m']
        supermarket_distances = df[df['kind'] == 'supermarket']['distance_to_duomo_m']
        control_distances = df[df['kind'] == 'control']['distance_to_duomo_m']
        
        if len(poi_distances) > 0:
            info(f"POI distances to Duomo (m): mean={poi_distances.mean():.0f}, min={poi_distances.min():.0f}, max={poi_distances.max():.0f}")
        
        if len(supermarket_distances) > 0:
            info(f"Supermarket distances to Duomo (m): mean={supermarket_distances.mean():.0f}, min={supermarket_distances.min():.0f}, max={supermarket_distances.max():.0f}")
        
        if len(control_distances) > 0:
            info(f"Control (POC) distances to Duomo (m): mean={control_distances.mean():.0f}, min={control_distances.min():.0f}, max={control_distances.max():.0f}")
        
        # Timestamps available
        timestamp_available = df['timestamp'].notna().sum()
        info(f"Features with timestamps: {timestamp_available}/{len(df)} ({timestamp_available/len(df)*100:.1f}%)")
        
        if timestamp_available > 0:
            df_with_ts = df[df['timestamp'].notna()].copy()
            df_with_ts['year'] = pd.to_datetime(df_with_ts['timestamp']).dt.year
            year_range = df_with_ts['year'].agg(['min', 'max'])
            info(f"Timestamp range: {year_range['min']} - {year_range['max']}")

    banner("All Done!")
    success(f"DDD Analysis outputs saved to: {output_folder}")

if __name__ == "__main__":
    try:
        main()
    finally:
        # Ensure clean return to the calling shell
        try:
            sys.stdout.flush()
            sys.stderr.flush()
        except Exception:
            pass
        os._exit(0)
