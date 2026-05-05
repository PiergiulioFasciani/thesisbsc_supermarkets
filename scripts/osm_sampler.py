#!/usr/bin/env python3
"""
OSM History Sampler
Extracts establishments from OSM history file and tracks their lifecycle
Outputs: Node ID, Type, Name, Location, Opening Date, Closing Date
"""

import osmium
import pandas as pd
import geopandas as gpd
import yaml
from datetime import datetime, timezone
from pathlib import Path
import logging
import argparse
from collections import defaultdict
import hashlib
import math
import os
import urllib.request
from shapely.geometry import Point, Polygon
from shapely.ops import transform
import pyproj
import sys

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
    UNDERLINE = '\033[4m'

# Custom logging formatter with colors
class ColoredFormatter(logging.Formatter):
    FORMATS = {
        logging.DEBUG: Colors.CYAN + '%(levelname)s' + Colors.ENDC + ' - %(message)s',
        logging.INFO: Colors.GREEN + '✓' + Colors.ENDC + ' %(message)s',
        logging.WARNING: Colors.YELLOW + '⚠' + Colors.ENDC + ' %(message)s',
        logging.ERROR: Colors.RED + '✗' + Colors.ENDC + ' %(message)s',
        logging.CRITICAL: Colors.RED + Colors.BOLD + '✗✗' + Colors.ENDC + ' %(message)s'
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

# Setup logging with colors
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(ColoredFormatter())
logging.basicConfig(level=logging.INFO, handlers=[handler])
logger = logging.getLogger(__name__)


class EstablishmentExtractor(osmium.SimpleHandler):
    """Handler to extract establishments and track their lifecycle"""
    
    def __init__(self, categories, center_lat=None, center_lon=None, radius_meters=None,
                 start_date=None, end_date=None):
        super().__init__()
        self.categories = categories  # Dict with pois, supermarkets, pocs
        self.node_versions = defaultdict(list)  # node_id -> list of versions
        self.center_lat = center_lat
        self.center_lon = center_lon
        self.radius_meters = radius_meters
        self.start_date = start_date
        self.end_date = end_date
        self.processed_count = 0
        self.processed_count = 0
    
    def _calculate_distance(self, lat1, lon1, lat2, lon2):
        """Calculate distance between two points using Haversine formula (in meters)"""
        R = 6371000  # Earth radius in meters
        phi1 = math.radians(lat1)
        phi2 = math.radians(lat2)
        delta_phi = math.radians(lat2 - lat1)
        delta_lambda = math.radians(lon2 - lon1)
        
        a = math.sin(delta_phi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda/2)**2
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        
        return R * c
    
    def _in_aoi(self, lat, lon):
        """Check if point is within Area of Interest (circular)"""
        if self.center_lat and self.center_lon and self.radius_meters:
            distance = self._calculate_distance(self.center_lat, self.center_lon, lat, lon)
            return distance <= self.radius_meters
        return True
    
    def _matches_category(self, tags):
        """Check if node matches any category and return (type, category_name)"""
        # Check POIs
        for category_name, tag_list in self.categories.get('pois', {}).items():
            for tag_filter in tag_list:
                if all(tags.get(k) == v for k, v in tag_filter.items()):
                    return ('POI', category_name)
        
        # Check Supermarkets
        for category_name, tag_list in self.categories.get('supermarkets', {}).items():
            for tag_filter in tag_list:
                if all(tags.get(k) == v for k, v in tag_filter.items()):
                    return ('Supermarket', category_name)
        
        # Check POCs
        for category_name, tag_list in self.categories.get('pocs', {}).items():
            for tag_filter in tag_list:
                if all(tags.get(k) == v for k, v in tag_filter.items()):
                    return ('POC', category_name)
        
        return None
        
    def node(self, n):
        """Process each node version"""
        try:
            # Filter by location
            if not self._in_aoi(n.location.lat, n.location.lon):
                return
            
            # Filter by date range
            if self.start_date and n.timestamp < self.start_date:
                return
            if self.end_date and n.timestamp > self.end_date:
                return
            
            # Extract tags
            tags = {tag.k: tag.v for tag in n.tags}
            
            # Check if node matches our categories
            match = self._matches_category(tags)
            if not match:
                return
            
            establishment_type, category_name = match
            
            # Store version data (minimal to save memory)
            version_data = {
                'node_id': n.id,
                'version': n.version,
                'timestamp': n.timestamp,
                'lon': n.location.lon,
                'lat': n.location.lat,
                'visible': n.visible,
                'type': establishment_type,
                'category': category_name,
                'name': tags.get('name', tags.get('brand', ''))
            }
            
            self.node_versions[n.id].append(version_data)
            self.processed_count += 1
            
            if self.processed_count % 10000 == 0:
                sys.stdout.write(f"\r{Colors.CYAN}⟳{Colors.ENDC} Processing: {self.processed_count:,} matching nodes found...")
                sys.stdout.flush()
            
        except Exception as e:
            logger.error(f"Error processing node {n.id}: {e}")
        except Exception as e:
            logger.error(f"Error processing node {n.id}: {e}")
    
    def _generate_node_id(self, osm_id, establishment_type):
        """Generate alphanumeric ID for establishment"""
        # Create a hash-based ID with type prefix
        type_prefix = {'POI': 'P', 'POC': 'C', 'Supermarket': 'S'}
        prefix = type_prefix.get(establishment_type, 'X')
        
        # Use last 7 digits of OSM ID plus a hash
        hash_val = hashlib.md5(str(osm_id).encode()).hexdigest()[:4].upper()
        return f"{prefix}{osm_id % 10000000:07d}{hash_val}"
    
    def _detect_repurposing(self, versions):
        """Detect if establishment was repurposed based on name/category changes"""
        lifecycles = []
        current_lifecycle = None
        
        for version in versions:
            if not current_lifecycle:
                # Start first lifecycle
                current_lifecycle = {
                    'opening_date': version['timestamp'],
                    'closing_date': None,
                    'type': version['type'],
                    'category': version['category'],
                    'name': version['name'],
                    'lon': version['lon'],
                    'lat': version['lat'],
                    'visible': version['visible']
                }
            else:
                # Check for major changes (simplified for speed)
                name_changed = (version['name'] != current_lifecycle['name'] and 
                               version['name'] and current_lifecycle['name'] and
                               len(version['name']) > 3 and len(current_lifecycle['name']) > 3)
                category_changed = version['category'] != current_lifecycle['category']
                
                if name_changed or category_changed:
                    # Close current lifecycle
                    current_lifecycle['closing_date'] = version['timestamp']
                    lifecycles.append(current_lifecycle)
                    
                    # Start new lifecycle
                    current_lifecycle = {
                        'opening_date': version['timestamp'],
                        'closing_date': None,
                        'type': version['type'],
                        'category': version['category'],
                        'name': version['name'],
                        'lon': version['lon'],
                        'lat': version['lat'],
                        'visible': version['visible']
                    }
                else:
                    # Update current lifecycle
                    current_lifecycle['closing_date'] = version['timestamp']
                    current_lifecycle['visible'] = version['visible']
                    current_lifecycle['lon'] = version['lon']
                    current_lifecycle['lat'] = version['lat']
        
        # Close final lifecycle
        if current_lifecycle:
            lifecycles.append(current_lifecycle)
        
        return lifecycles
    
    def get_dataframe(self):
        """Convert collected nodes to DataFrame with lifecycle tracking"""
        if not self.node_versions:
            logger.warning("No establishments extracted")
            return pd.DataFrame()
        
        logger.info(f"Processing {len(self.node_versions)} unique nodes...")
        
        establishments = []
        
        for osm_id, versions in self.node_versions.items():
            # Sort versions by timestamp
            versions.sort(key=lambda x: x['timestamp'])
            
            # Detect lifecycles (opening/closing/repurposing)
            lifecycles = self._detect_repurposing(versions)
            
            # Create entry for each lifecycle
            for i, lifecycle in enumerate(lifecycles):
                node_id = self._generate_node_id(osm_id, lifecycle['type'])
                if i > 0:
                    node_id += f"R{i}"  # Add repurposing suffix
                
                establishment = {
                    'node_id': node_id,
                    'osm_id': osm_id,
                    'type': lifecycle['type'],
                    'category': lifecycle['category'],
                    'name': lifecycle['name'] if lifecycle['name'] else 'Unknown',
                    'latitude': lifecycle['lat'],
                    'longitude': lifecycle['lon'],
                    'opening_date': lifecycle['opening_date'].strftime('%Y-%m-%d'),
                    'closing_date': lifecycle['closing_date'].strftime('%Y-%m-%d') if lifecycle['closing_date'] else None,
                    'still_open': lifecycle['visible'] and lifecycle['closing_date'] is None
                }
                
                establishments.append(establishment)
        
        df = pd.DataFrame(establishments)
        
        # Sort by type, then opening date
        df = df.sort_values(['type', 'opening_date', 'node_id'])
        
        logger.info(f"Extracted {len(df)} establishment lifecycles from {len(self.node_versions)} unique nodes")
        
        return df


def load_config(config_file):
    """Load configuration from YAML file"""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config


def get_workflow_config(config, section_name):
    """Return a workflow-specific section or fall back to the root config."""
    section = config.get(section_name)
    if isinstance(section, dict):
        return section
    return config


def ensure_osm_file_exists(input_file, google_drive_file_id=None):
    """Ensure the OSM file exists locally, downloading it from Google Drive if needed."""
    input_path = Path(input_file)

    if input_path.exists():
        logger.info(f"OSM file found: {input_path}")
        return str(input_path)

    drive_id = google_drive_file_id or os.environ.get('OSM_GOOGLE_DRIVE_FILE_ID')
    if not drive_id:
        raise FileNotFoundError(
            f"OSM file not found: {input_path}\n"
            f"Set `sampler.google_drive_file_id` in config.yml or export OSM_GOOGLE_DRIVE_FILE_ID."
        )

    input_path.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://drive.google.com/uc?id={drive_id}&export=download&confirm=t"

    def reporthook(block_num, block_size, total_size):
        downloaded = block_num * block_size
        if total_size > 0:
            percent = min(downloaded * 100 / total_size, 100)
            mb_downloaded = downloaded / (1024 * 1024)
            mb_total = total_size / (1024 * 1024)
            print(f"\r  Download progress: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)", end='', flush=True)

    logger.info(f"Downloading OSM file from Google Drive (ID: {drive_id[:10]}...)")
    urllib.request.urlretrieve(url, str(input_path), reporthook=reporthook)
    print()

    if not input_path.exists():
        raise RuntimeError("Download completed but the file was not created")

    logger.info(f"Download complete: {input_path}")
    return str(input_path)


def create_aoi_polygon(center_lat, center_lon, radius_meters):
    """Create a circular polygon for the AOI"""
    # Create a circle using azimuthal equidistant projection
    local_azimuthal_projection = f"+proj=aeqd +R=6371000 +units=m +lat_0={center_lat} +lon_0={center_lon}"
    wgs84 = pyproj.CRS('EPSG:4326')
    aeqd = pyproj.CRS(local_azimuthal_projection)
    
    project = pyproj.Transformer.from_crs(wgs84, aeqd, always_xy=True).transform
    project_back = pyproj.Transformer.from_crs(aeqd, wgs84, always_xy=True).transform
    
    # Create circle in projected coordinates
    center_point = Point(0, 0)
    circle = center_point.buffer(radius_meters)
    
    # Transform back to WGS84
    circle_wgs84 = transform(project_back, circle)
    
    return circle_wgs84


def extract_establishments(input_file, output_dir, config_file=None, 
                          center_lat=None, center_lon=None, radius_meters=None,
                          start_date=None, end_date=None):
    """
    Extract establishments from OSM history file with lifecycle tracking
    
    Args:
        input_file: Path to .osh.pbf file
        output_dir: Directory for output files (will create timestamped subfolder)
        config_file: Optional YAML config file
        center_lat, center_lon, radius_meters: Circular AOI
        start_date: Start date as ISO string
        end_date: End date as ISO string
    """
    # Load config if provided
    categories = {'pois': {}, 'supermarkets': {}, 'pocs': {}}
    google_drive_file_id = None
    
    if config_file:
        config = load_config(config_file)
        sampler_config = get_workflow_config(config, 'sampler')
        categories = {
            'pois': sampler_config.get('pois', {}),
            'supermarkets': sampler_config.get('supermarkets', {}),
            'pocs': sampler_config.get('pocs', {})
        }
        google_drive_file_id = sampler_config.get('google_drive_file_id')
        
        # Get AOI parameters from config
        if not center_lat:
            center_lat = sampler_config.get('center_lat')
            center_lon = sampler_config.get('center_lon')
            radius_meters = sampler_config.get('radius_meters')
        
        if not start_date:
            start_date = sampler_config.get('start_date')
        if not end_date:
            end_date = sampler_config.get('end_date')

    input_file = ensure_osm_file_exists(input_file, google_drive_file_id)
    
    # Parse dates to timezone-aware datetime
    start_date_parsed = None
    end_date_parsed = None
    if start_date:
        start_date_parsed = datetime.fromisoformat(start_date).replace(tzinfo=timezone.utc)
    if end_date:
        end_date_parsed = datetime.fromisoformat(end_date).replace(tzinfo=timezone.utc)
    # Create timestamped output folder
    timestamp = datetime.now().strftime('%d%m%y_%H%M%S')
    sample_dir = Path(output_dir) / 'samples' / f'sample_{timestamp}'
    sample_dir.mkdir(parents=True, exist_ok=True)
    
    # Print beautiful header
    print("\n" + "="*80)
    print(f"{Colors.BOLD}{Colors.HEADER}OSM ESTABLISHMENT SAMPLER{Colors.ENDC}")
    print("="*80)
    print(f"\n{Colors.BOLD}Configuration:{Colors.ENDC}")
    print(f"  Input file: {Colors.CYAN}{Path(input_file).name}{Colors.ENDC}")
    print(f"  AOI Center: {Colors.CYAN}({center_lat:.4f}, {center_lon:.4f}){Colors.ENDC}")
    print(f"  Radius: {Colors.CYAN}{radius_meters:,}m{Colors.ENDC} ({radius_meters/1000:.1f} km)")
    if start_date or end_date:
        print(f"  Date range: {Colors.CYAN}{start_date or 'any'}{Colors.ENDC} to {Colors.CYAN}{end_date or 'any'}{Colors.ENDC}")
    print(f"  Output: {Colors.CYAN}{sample_dir.name}/{Colors.ENDC}")
    print("\n" + "-"*80 + "\n")
    
    # Create handler
    handler = EstablishmentExtractor(
        categories=categories,
        center_lat=center_lat,
        center_lon=center_lon,
        radius_meters=radius_meters,
        start_date=start_date_parsed,
        end_date=end_date_parsed
    )
    
    # Process the file
    print(f"{Colors.BOLD}[1/4] Extracting from OSM history...{Colors.ENDC}")
    handler.apply_file(input_file)
    print(f"\r{Colors.GREEN}✓{Colors.ENDC} Extracted {handler.processed_count:,} matching node versions" + " "*20)
    
    # Get results as DataFrame
    print(f"\n{Colors.BOLD}[2/4] Processing lifecycles...{Colors.ENDC}")
    df = handler.get_dataframe()
    
    if df.empty:
        print(f"\n{Colors.RED}✗ No data extracted{Colors.ENDC}")
        return
    
    # Save CSV
    print(f"\n{Colors.BOLD}[3/4] Saving data...{Colors.ENDC}")
    csv_file = sample_dir / 'establishments.csv'
    df.to_csv(csv_file, index=False)
    logger.info(f"CSV: {len(df):,} records → {csv_file.name}")
    
    # Create GeoDataFrame for spatial output
    gdf = gpd.GeoDataFrame(
        df,
        geometry=[Point(lon, lat) for lon, lat in zip(df['longitude'], df['latitude'])],
        crs='EPSG:4326'
    )
    
    # Create AOI polygon
    aoi_polygon = create_aoi_polygon(center_lat, center_lon, radius_meters)
    aoi_gdf = gpd.GeoDataFrame(
        [{'name': 'AOI', 'center_lat': center_lat, 'center_lon': center_lon, 
          'radius_m': radius_meters}],
        geometry=[aoi_polygon],
        crs='EPSG:4326'
    )
    
    # Save to GeoPackage with multiple layers
    print(f"{Colors.BOLD}[4/4] Creating GeoPackage layers...{Colors.ENDC}")
    gpkg_file = sample_dir / 'establishments.gpkg'
    aoi_gdf.to_file(gpkg_file, layer='aoi', driver='GPKG')
    gdf.to_file(gpkg_file, layer='establishments', driver='GPKG')
    
    layer_count = 2  # aoi + establishments
    # Create only open/closed layers (skip individual type layers for speed)
    open_data = gdf[gdf['still_open'] == True]
    closed_data = gdf[gdf['still_open'] == False]
    
    if not open_data.empty:
        open_data.to_file(gpkg_file, layer='open', driver='GPKG')
        layer_count += 1
    
    if not closed_data.empty:
        closed_data.to_file(gpkg_file, layer='closed', driver='GPKG')
        layer_count += 1
    
    logger.info(f"GeoPackage: {layer_count} layers → {gpkg_file.name}")
    
    # Create metadata file
    metadata_file = sample_dir / 'metadata.txt'
    with open(metadata_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("OSM ESTABLISHMENT EXTRACTION METADATA\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Extraction Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input File: {input_file}\n\n")
        
        f.write("AREA OF INTEREST (AOI):\n")
        f.write(f"  Center Latitude: {center_lat}\n")
        f.write(f"  Center Longitude: {center_lon}\n")
        f.write(f"  Radius: {radius_meters} meters ({radius_meters/1000:.2f} km)\n\n")
        
        f.write("DATE RANGE:\n")
        f.write(f"  Start Date: {start_date if start_date else 'Not specified'}\n")
        f.write(f"  End Date: {end_date if end_date else 'Not specified'}\n\n")
        
        f.write("EXTRACTION RESULTS:\n")
        f.write(f"  Total Establishments: {len(df)}\n")
        f.write(f"  Unique OSM Nodes: {df['osm_id'].nunique()}\n\n")
        
        f.write("BY TYPE:\n")
        for est_type in ['POI', 'Supermarket', 'POC']:
            count = len(df[df['type'] == est_type])
            open_count = len(df[(df['type'] == est_type) & (df['still_open'] == True)])
            closed_count = count - open_count
            if count > 0:
                f.write(f"  {est_type}:\n")
                f.write(f"    Total: {count}\n")
                f.write(f"    Open: {open_count}\n")
                f.write(f"    Closed: {closed_count}\n")
        
        f.write(f"\nDATE RANGE IN DATA:\n")
        f.write(f"  Earliest Opening: {df['opening_date'].min()}\n")
        f.write(f"  Latest Opening: {df['opening_date'].max()}\n")
        f.write(f"  Still Open: {df['still_open'].sum()}\n")
        f.write(f"  Closed: {(~df['still_open']).sum()}\n\n")
        
        f.write("TOP CATEGORIES:\n")
        for category, count in df['category'].value_counts().head(15).items():
            f.write(f"  {category}: {count}\n")
        
        f.write(f"\nCATEGORY DEFINITIONS:\n")
        f.write("\nPOIs (Points of Interest):\n")
        for cat, tags in categories.get('pois', {}).items():
            f.write(f"  {cat}: {tags}\n")
        
        f.write("\nSupermarkets:\n")
        for cat, tags in categories.get('supermarkets', {}).items():
            f.write(f"  {cat}: {tags}\n")
        
        f.write("\nPOCs (Points of Comparison):\n")
    logger.info(f"Metadata: {metadata_file.name}")
    
    # Print beautiful summary
    print("\n" + "="*80)
    print(f"{Colors.BOLD}{Colors.GREEN}✓ EXTRACTION COMPLETE{Colors.ENDC}")
    print("="*80)
    
    print(f"\n{Colors.BOLD}Results:{Colors.ENDC}")
    print(f"  Total establishments: {Colors.CYAN}{len(df):,}{Colors.ENDC}")
    print(f"  Unique OSM nodes: {Colors.CYAN}{df['osm_id'].nunique():,}{Colors.ENDC}")
    print(f"  Still open: {Colors.GREEN}{df['still_open'].sum():,}{Colors.ENDC}")
    print(f"  Closed: {Colors.YELLOW}{(~df['still_open']).sum():,}{Colors.ENDC}")
    
    print(f"\n{Colors.BOLD}By Type:{Colors.ENDC}")
    for est_type in ['POI', 'Supermarket', 'POC']:
        count = len(df[df['type'] == est_type])
        if count > 0:
            open_count = len(df[(df['type'] == est_type) & (df['still_open'] == True)])
            closed_count = count - open_count
            print(f"  {est_type:12} {Colors.CYAN}{count:5,}{Colors.ENDC} total  " +
                  f"({Colors.GREEN}{open_count:,}{Colors.ENDC} open, " +
                  f"{Colors.YELLOW}{closed_count:,}{Colors.ENDC} closed)")
    
    print(f"\n{Colors.BOLD}Top Categories:{Colors.ENDC}")
    for category, count in df['category'].value_counts().head(5).items():
        print(f"  {category:25} {Colors.CYAN}{count:5,}{Colors.ENDC}")
    
    print(f"\n{Colors.BOLD}Output:{Colors.ENDC} {Colors.CYAN}{sample_dir}{Colors.ENDC}")
    print(f"  {Colors.GREEN}→{Colors.ENDC} establishments.csv")
    print(f"  {Colors.GREEN}→{Colors.ENDC} establishments.gpkg ({layer_count} layers)")
    print(f"  {Colors.GREEN}→{Colors.ENDC} metadata.txt")
    print("="*80 + "\n")
    # Print summary to console
    print("\n" + "="*80)
    print("EXTRACTION SUMMARY")
    print("="*80)
    print(f"Total establishments: {len(df)}")
    print(f"Unique OSM nodes: {df['osm_id'].nunique()}")
    print(f"\nBy Type:")
    for est_type in ['POI', 'Supermarket', 'POC']:
        count = len(df[df['type'] == est_type])
        if count > 0:
            print(f"  {est_type}: {count}")
    print(f"\nOutput directory: {sample_dir}")
    print(f"  - establishments.csv")
    print(f"  - establishments.gpkg")
    print(f"  - metadata.txt")
    print("="*80)


def main():
    parser = argparse.ArgumentParser(
        description='Extract establishments from OSM history file with lifecycle tracking',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using config file
    python osm_sampler.py --config config.yml
  
  # Direct parameters with circular AOI
  python osm_sampler.py --input history.osh.pbf --output-dir data/output \\
      --center-lat 45.4642 --center-lon 9.1900 --radius 12000 \\
      --start-date 2015-01-01 --end-date 2024-12-31

Config file format (YAML):
    sampler:
        input_file: data/input/history.osh.pbf
        output_dir: data/output
        center_lat: 45.4642
        center_lon: 9.1900
        radius_meters: 12000
        start_date: "2015-01-01"
        end_date: "2024-12-31"
        pois:
            bakeries:
                - shop: bakery
        supermarkets:
            supermarkets:
                - shop: supermarket
        pocs:
            restaurants_fast_food:
                - amenity: restaurant
        """
    )
    
    parser.add_argument('--config', '-c', type=str,
                       help='YAML configuration file')
    parser.add_argument('--input', '-i', type=str,
                       help='Input .osh.pbf file')
    parser.add_argument('--output-dir', '-o', type=str,
                       help='Output directory (will create timestamped subfolder)')
    parser.add_argument('--center-lat', type=float,
                       help='Center latitude for circular AOI')
    parser.add_argument('--center-lon', type=float,
                       help='Center longitude for circular AOI')
    parser.add_argument('--radius', type=float,
                       help='Radius in meters for circular AOI')
    parser.add_argument('--start-date', type=str,
                       help='Start date (ISO format: YYYY-MM-DD)')
    parser.add_argument('--end-date', type=str,
                       help='End date (ISO format: YYYY-MM-DD)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.config:
        config = load_config(args.config)
        sampler_config = get_workflow_config(config, 'sampler')
        input_file = sampler_config.get('input_file')
        output_dir = sampler_config.get('output_dir', 'data/output')
        extract_establishments(input_file, output_dir, config_file=args.config)
    elif args.input and args.output_dir:
        extract_establishments(
            input_file=args.input,
            output_dir=args.output_dir,
            center_lat=args.center_lat,
            center_lon=args.center_lon,
            radius_meters=args.radius,
            start_date=args.start_date,
            end_date=args.end_date
        )
    else:
        parser.error("Either --config or both --input and --output-dir are required")


if __name__ == '__main__':
    main()
