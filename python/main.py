# -------------------------------------------------------------------------------------------------
# Quadtree tiles → ONE CSV (lean, interpretable headers)
#
# Tile-level distances (meters), with 0→outside fallback:
#   d_pt_m     : tile polygon → nearest public-transport stop distance (m);
#                if 0 (a stop lies inside the tile), uses nearest OUTSIDE-OF-TILE stop.
#   d_super_m  : tile polygon → nearest supermarket distance (m);
#                if 0 (a supermarket lies inside the tile), uses nearest OUTSIDE-OF-TILE supermarket.
#   d_duomo_m  : tile polygon → Duomo distance (m);
#                if Duomo lies inside the tile (distance==0), uses tile-centroid → Duomo.
#
# Buffer counts (polygon-based; unchanged semantics, compact names):
#   pt200_cnt, pt400_cnt, super300_cnt, super600_cnt
#
# Counts inside tile (unchanged):
#   n_pois, n_poc + has_poi/has_poc
#
# CLI (the run_all.sh/.bat launchers pass these)
#   python python/main.py --center "45.43318,9.18378" --radius-m 1200 \
#       --pbf data/pbf/nord-ovest-latest.osm.pbf --out-dir data --min-tile-m 50
#
# Requirements
#   pip install pyosmium geopandas shapely pyproj fiona pandas numpy tqdm rich
# -------------------------------------------------------------------------------------------------

from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon, box
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

# Target POIs (drive the quadtree)
POI_SHOPS = {"bakery", "greengrocer", "butcher", "convenience"}

# Supermarkets
SUPERMARKET_SHOPS = {"supermarket"}

# Public transport
PT_HIGHWAY = {"bus_stop"}
PT_AMENITY = {"bus_station"}
PT_RAILWAY = {"station", "halt", "tram_stop", "subway_entrance", "stop", "stop_position", "platform"}
PT_PUBLIC_TRANSPORT = {"station", "platform", "stop_position", "stop_area"}

# Control (POC): mortuary-related
MORTUARY_SHOP = {"funeral_directors", "funeral_home"}
MORTUARY_AMENITY = {"crematorium", "grave_yard"}
MORTUARY_LANDUSE = {"cemetery"}

# ---- Portable defaults for run_all.sh/.bat (cross-platform) ----
DEFAULT_PBF = os.path.join("data", "pbf", "nord-ovest-latest.osm.pbf")
DEFAULT_OUT_BASE = "data"

# Buffer radii (meters) for counts
PT_BUFFER_LIST = [200, 400]
SUPER_BUFFER_LIST = [300, 600]

# ----------------------------------------- UX -----------------------------------------

def banner(title: str):
    if RICH: console.print(Panel.fit(f"[bold cyan]{title}[/]", border_style="cyan"))
    else: print(f"\n=== {title} ===")

def info(msg: str):
    if RICH: console.print(f"[bold]• {msg}[/]")
    else: print(f"• {msg}")

def success(msg: str):
    if RICH: console.print(f"[green]✓ {msg}[/]")
    else: print(f"[OK] {msg}")

def warn(msg: str):
    if RICH: console.print(f"[yellow]! {msg}[/]")
    else: print(f"[WARN] {msg}")

def err(msg: str):
    if RICH: console.print(f"[red]✗ {msg}[/]")
    else: print(f"[ERR] {msg}")

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
    to_wgs  = Transformer.from_crs(aeqd, wgs84, always_xy=True).transform
    poly_m = shp_transform(to_aeqd, Point(center_lon, center_lat)).buffer(radius_m, resolution=64)
    poly   = shp_transform(to_wgs, poly_m)
    try: poly = poly.buffer(0)
    except Exception: pass
    return poly

def local_utm_epsg(lon: float, lat: float) -> int:
    zone = int((lon + 180.0) / 6.0) + 1
    return (32600 if lat >= 0 else 32700) + zone

def fiona_safe(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    if gdf is None or len(gdf) == 0: return gdf
    gdf = gdf.copy()
    for col in gdf.columns:
        if col == gdf.geometry.name: continue
        dt = gdf[col].dtype
        if "Int" in str(dt): gdf[col] = gdf[col].astype(np.int64)
        elif pd.api.types.is_string_dtype(dt): gdf[col] = gdf[col].astype(object)
        elif pd.api.types.is_bool_dtype(dt): gdf[col] = gdf[col].astype(np.int8)
    return gdf

def prepare_output_dir(base_dir: str, radius_m: float):
    """
    Create an output directory like:
      quadtree_YYYYMMDD_HHMMSS_<AOI meters>
    Example: quadtree_20250918_143918_1200
    """
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    radius_tag = f"{int(round(radius_m))}"
    out_dir = os.path.join(base_dir, f"quadtree_{ts}_{radius_tag}")
    os.makedirs(out_dir, exist_ok=True)
    return out_dir

# -------------------------- DISTANCES / COUNTS (poly only + fallback) --------------------------

def _points_inside_mask(poly_m: Polygon, points_gdf_m: gpd.GeoDataFrame) -> np.ndarray:
    if points_gdf_m is None or len(points_gdf_m) == 0:
        return np.zeros(0, dtype=bool)
    poly_series = gpd.GeoSeries([poly_m], crs=points_gdf_m.crs)
    return points_gdf_m.within(poly_series.iloc[0]).to_numpy()

def min_poly_to_nearest_point_any(poly_m, points_gdf_m, sidx, steps=(50,100,200,400,800,1600,3200,6400)) -> float:
    if points_gdf_m is None or len(points_gdf_m) == 0:
        return float("nan")
    inside_mask = _points_inside_mask(poly_m, points_gdf_m)
    if inside_mask.size and inside_mask.any():
        return 0.0
    for r in steps:
        buf = poly_m.buffer(r)
        cand_idx = list(sidx.query(buf, predicate="intersects"))
        if not cand_idx: continue
        pts = points_gdf_m.geometry.iloc[cand_idx]
        return float(min(poly_m.distance(p) for p in pts))
    return float(min(poly_m.distance(p) for p in points_gdf_m.geometry))

def min_poly_to_nearest_point_outside(poly_m, points_gdf_m, sidx, steps=(50,100,200,400,800,1600,3200,6400)) -> float:
    if points_gdf_m is None or len(points_gdf_m) == 0:
        return float("nan")
    inside_mask = _points_inside_mask(poly_m, points_gdf_m)
    if inside_mask.size == 0:
        return float("nan")
    if (~inside_mask).sum() == 0:
        return float("nan")
    for r in steps:
        buf = poly_m.buffer(r)
        cand_idx = list(sidx.query(buf, predicate="intersects"))
        if not cand_idx: continue
        outside_idx = [i for i in cand_idx if not inside_mask[i]]
        if not outside_idx: continue
        pts = points_gdf_m.geometry.iloc[outside_idx]
        return float(min(poly_m.distance(p) for p in pts))
    pts_out = points_gdf_m.geometry[~inside_mask]
    if len(pts_out) == 0: return float("nan")
    return float(min(poly_m.distance(p) for p in pts_out))

def count_points_within_buffer(poly_m, points_gdf_m: gpd.GeoDataFrame, sidx, radius_m: float) -> int:
    if points_gdf_m is None or len(points_gdf_m) == 0:
        return 0
    buf = poly_m.buffer(radius_m)
    cand = list(sidx.query(buf, predicate="intersects"))
    if not cand:
        return 0
    sub = points_gdf_m.iloc[cand]
    return int(sub.within(gpd.GeoSeries([buf], crs=points_gdf_m.crs).iloc[0]).sum())

# ------------------------------------- OSM PARSER -------------------------------------

@dataclass
class Feature:
    id: int
    lon: float
    lat: float
    name: Optional[str]
    kind: str          # 'poi_*', 'supermarket', 'transport', 'control_mortuary'
    source: str        # 'node' | 'area'

class TagMatcher:
    def match_kind(self, tags) -> Optional[str]:
        shop = tags.get("shop")
        amenity = tags.get("amenity")
        railway = tags.get("railway")
        pubtr = tags.get("public_transport")
        landuse = tags.get("landuse")

        if shop in POI_SHOPS: return f"poi_{shop}"
        if shop in SUPERMARKET_SHOPS: return "supermarket"
        if (shop in MORTUARY_SHOP) or (amenity in MORTUARY_AMENITY) or (landuse in MORTUARY_LANDUSE):
            return "control_mortuary"
        if (tags.get("highway") in PT_HIGHWAY or
            amenity in PT_AMENITY or
            railway in PT_RAILWAY or
            pubtr in PT_PUBLIC_TRANSPORT):
            return "transport"
        return None

class AOIFilter:
    def __init__(self, aoi_poly: Polygon):
        self.aoi = aoi_poly
        self.bounds = aoi_poly.bounds
    def contains_lonlat(self, lon: float, lat: float) -> bool:
        minx, miny, maxx, maxy = self.bounds
        if not (minx <= lon <= maxx and miny <= lat <= maxy):
            return False
        return self.aoi.covers(Point(lon, lat))

class EnhancedAOIFilter:
    """Enhanced AOI filter with expanded bounding box for global collections"""
    def __init__(self, aoi_poly: Polygon, minx: float, miny: float, maxx: float, maxy: float):
        self.aoi = aoi_poly
        self.aoi_bounds = aoi_poly.bounds
        self.expanded_bounds = (minx, miny, maxx, maxy)
        
    def contains_lonlat(self, lon: float, lat: float) -> bool:
        """Check if point is in AOI (for AOI-only features)"""
        minx, miny, maxx, maxy = self.aoi_bounds
        if not (minx <= lon <= maxx and miny <= lat <= maxy):
            return False
        return self.aoi.covers(Point(lon, lat))
    
    def in_expanded_area(self, lon: float, lat: float) -> bool:
        """Check if point is in expanded area (for global collections)"""
        minx, miny, maxx, maxy = self.expanded_bounds
        return (minx <= lon <= maxx and miny <= lat <= maxy)

class OSMCollector(osm.SimpleHandler):
    """ AOI features: POIs, control (mortuary), supermarkets & PT (for GIS);
        GLOBAL sets: supermarkets + PT (for distances & buffer counts) """
    def __init__(self, tm: TagMatcher, af):
        super().__init__()
        self.tm, self.af = tm, af
        self.wkbf = osm.geom.WKBFactory()
        self.features: List[Feature] = []
        self.supers_global: List[Feature] = []
        self.pt_global: List[Feature] = []

    def _maybe_keep_aoi(self, f: Feature):
        if self.af.contains_lonlat(f.lon, f.lat):
            self.features.append(f)
    
    def _maybe_keep_global(self, f: Feature):
        """Keep for global collections if in expanded area"""
        if hasattr(self.af, 'in_expanded_area'):
            return self.af.in_expanded_area(f.lon, f.lat)
        else:
            return True  # fallback for regular AOIFilter

    def node(self, n):
        kind = self.tm.match_kind(n.tags)
        if kind is None or not n.location.valid(): return
        lon, lat = n.location.lon, n.location.lat
        
        # Skip if completely outside expanded area for efficiency
        if hasattr(self.af, 'in_expanded_area') and not self.af.in_expanded_area(lon, lat):
            return
            
        name = n.tags.get("name")
        if kind == "supermarket":
            if self._maybe_keep_global(Feature(n.id, lon, lat, name, kind, "node")):
                self.supers_global.append(Feature(n.id, lon, lat, name, kind, "node"))
            self._maybe_keep_aoi(Feature(n.id, lon, lat, name, kind, "node")); return
        if kind == "transport":
            if self._maybe_keep_global(Feature(n.id, lon, lat, name, kind, "node")):
                self.pt_global.append(Feature(n.id, lon, lat, name, kind, "node"))
            self._maybe_keep_aoi(Feature(n.id, lon, lat, name, "transport", "node")); return
        self._maybe_keep_aoi(Feature(n.id, lon, lat, name, kind, "node"))

    def area(self, a):
        kind = self.tm.match_kind(a.tags)
        if kind is None: return
        try:
            wkb = None
            try: wkb = self.wkbf.create_multipolygon(a)
            except Exception: pass
            if wkb is None:
                try: wkb = self.wkbf.create_polygon(a)
                except Exception: return
            geom = shp_wkb.loads(wkb, hex=True)
            c = geom.centroid
            lon, lat = c.x, c.y
            
            # Skip if completely outside expanded area for efficiency
            if hasattr(self.af, 'in_expanded_area') and not self.af.in_expanded_area(lon, lat):
                return
                
            name = a.tags.get("name")
            if kind == "supermarket":
                if self._maybe_keep_global(Feature(a.id, lon, lat, name, kind, "area")):
                    self.supers_global.append(Feature(a.id, lon, lat, name, kind, "area"))
                self._maybe_keep_aoi(Feature(a.id, lon, lat, name, kind, "area")); return
            if kind == "transport":
                if self._maybe_keep_global(Feature(a.id, lon, lat, name, kind, "area")):
                    self.pt_global.append(Feature(a.id, lon, lat, name, kind, "area"))
                self._maybe_keep_aoi(Feature(a.id, lon, lat, name, "transport", "area")); return
            self._maybe_keep_aoi(Feature(a.id, lon, lat, name, kind, "area"))
        except Exception:
            return

def parse_pbf_inside_aoi(pbf_path: str, aoi_poly: Polygon):
    banner("Step 2 — Reading OSM (AOI parse + GLOBAL supers/PT with efficient filtering)")
    tm, af = TagMatcher(), AOIFilter(aoi_poly)
    
    # Get bounding box for efficient pre-filtering
    minx, miny, maxx, maxy = aoi_poly.bounds
    # Expand by ~2km buffer to ensure we get all global supermarkets/PT for distance calculations
    buffer_deg = 0.02  # approximately 2km at latitude ~45°
    parse_minx, parse_miny = minx - buffer_deg, miny - buffer_deg
    parse_maxx, parse_maxy = maxx + buffer_deg, maxy + buffer_deg
    info(f"Target area: lon [{parse_minx:.6f}, {parse_maxx:.6f}], lat [{parse_miny:.6f}, {parse_maxy:.6f}]")
    
    # Use enhanced filter that pre-screens by bounding box
    enhanced_af = EnhancedAOIFilter(aoi_poly, parse_minx, parse_miny, parse_maxx, parse_maxy)
    
    if RICH:
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task = progress.add_task("Parsing PBF…", total=None)
            h = OSMCollector(tm, enhanced_af)
            h.apply_file(str(pbf_path), locations=True, idx="flex_mem")
            progress.update(task, description="Finalizing…")
    else:
        h = OSMCollector(tm, enhanced_af)
        h.apply_file(str(pbf_path), locations=True, idx="flex_mem")

    def to_gdf(items: List[Feature]):
        if not items:
            return gpd.GeoDataFrame(columns=["kind","name","source","osm_id","geometry"], crs="EPSG:4326")
        df = pd.DataFrame([(f.lon, f.lat, f.kind, f.name, f.source, f.id) for f in items],
                          columns=["lon","lat","kind","name","source","osm_id"])
        return gpd.GeoDataFrame(df[["kind","name","source","osm_id"]],
                                geometry=[Point(lon, lat) for lon, lat, *_ in df.values],
                                crs="EPSG:4326")

    all_aoi = to_gdf(h.features)
    pois_wgs = all_aoi[all_aoi["kind"].str.startswith("poi_")].copy()
    poc_wgs  = all_aoi[all_aoi["kind"]=="control_mortuary"].copy()
    supers_wgs = all_aoi[all_aoi["kind"]=="supermarket"].copy()
    pt_wgs     = all_aoi[all_aoi["kind"]=="transport"].copy()

    supers_global_wgs = to_gdf(h.supers_global)
    pt_global_wgs     = to_gdf(h.pt_global)

    info(f"AOI kept: {len(all_aoi)} | POIs: {len(pois_wgs)}, POCs: {len(poc_wgs)}, Supers(AOI): {len(supers_wgs)}, PT(AOI): {len(pt_wgs)}")
    info(f"GLOBAL — Supermarkets: {len(supers_global_wgs)}, PT: {len(pt_global_wgs)}")
    return pois_wgs, poc_wgs, supers_wgs, pt_wgs, supers_global_wgs, pt_global_wgs

# ------------------------------------- QUADTREE -------------------------------------

def split_into_quadrants(poly):
    minx, miny, maxx, maxy = poly.bounds
    cx, cy = (minx + maxx)/2.0, (miny + maxy)/2.0
    return [box(minx, miny, cx, cy), box(cx, miny, maxx, cy),
            box(minx, cy, cx, maxy), box(cx, cy, maxx, maxy)]

def build_quadtree(root_square_m: Polygon, poi_points_m: gpd.GeoDataFrame, min_tile_m: float):
    banner("Step 4 — Building quadtree (split iff ≥ 2 target POIs)")
    if len(poi_points_m) == 0:
        warn("No target POIs in AOI; grid will be one tile.")
        return [(root_square_m, 0, 0)]

    sindex = poi_points_m.sindex
    side_root = root_square_m.bounds[2] - root_square_m.bounds[0]

    def count_in(poly_m: Polygon) -> int:
        idx = list(sindex.query(poly_m, predicate="intersects"))
        if not idx: return 0
        sub = poi_points_m.iloc[idx]
        return int(sub.within(gpd.GeoSeries([poly_m], crs=poi_points_m.crs).iloc[0]).sum())

    leaves, stack, visited = [], [root_square_m], 0
    if RICH:
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
                      TimeElapsedColumn()) as progress:
            task = progress.add_task("Splitting tiles…", total=None)
            while stack:
                tile = stack.pop()
                visited += 1
                side = tile.bounds[2] - tile.bounds[0]
                depth = 0 if side == 0 else int(round(math.log2(side_root / side)))
                n = count_in(tile)
                progress.update(task, description=f"Splitting tiles… visited={visited:,} leaves={len(leaves):,} depth={depth} n={n}")
                if (n >= 2) and (side > min_tile_m) and (len(leaves) + len(stack) < 2_000_000):
                    stack.extend(split_into_quadrants(tile))
                else:
                    leaves.append((tile, depth, n))
    else:
        while stack:
            tile = stack.pop()
            visited += 1
            side = tile.bounds[2] - tile.bounds[0]
            depth = 0 if side == 0 else int(round(math.log2(side_root / side)))
            n = count_in(tile)
            if (n >= 2) and (side > min_tile_m) and (len(leaves) + len(stack) < 2_000_000):
                stack.extend(split_into_quadrants(tile))
            else:
                leaves.append((tile, depth, n))
    success(f"Quadtree built: leaves={len(leaves):,}, visited={visited:,}")
    return leaves

# ---------------------------------------- MAIN ----------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Quadtree tiles with lean CSV (outside fallback for 0 distances)")
    ap.add_argument("--center", required=True, help="lat,lon (e.g. '45.43318,9.18378')")
    ap.add_argument("--radius-m", type=float, required=True)
    ap.add_argument("--pbf", default=DEFAULT_PBF)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_BASE)
    ap.add_argument("--min-tile-m", type=float, default=50.0)
    args = ap.parse_args()

    # Step 1 — Setup
    banner("Step 1 — Setup")
    center_lat, center_lon = parse_center(args.center)
    pbf_path = Path(args.pbf).resolve()
    if not pbf_path.exists():
        err(f"PBF not found: {pbf_path}")
        raise SystemExit(1)

    out_dir = prepare_output_dir(args.out_dir, args.radius_m)
    info(f"Center: ({center_lat:.6f}, {center_lon:.6f})  Radius: {args.radius_m:.0f} m")
    info(f"PBF   : {pbf_path}")
    info(f"Out   : {out_dir}")

    # Build AOI polygon first (WGS84)
    banner("Step 1a — Building AOI")
    aoi_poly = circle_polygon_wgs84(center_lon=center_lon, center_lat=center_lat, radius_m=args.radius_m)
    minx, miny, maxx, maxy = aoi_poly.bounds
    info(f"AOI bbox (WGS84): lon [{minx:.6f}, {maxx:.6f}], lat [{miny:.6f}, {maxy:.6f}]")
    
    # Calculate area reduction for efficiency info
    total_area_km2 = (maxx - minx) * (maxy - miny) * 111 * 111  # rough km²
    info(f"AOI area: ~{abs(total_area_km2):.1f} km² (parsing with ~2km buffer)")

    # Step 2 — Parse OSM
    pois_wgs, poc_wgs, supers_wgs, pt_wgs, supers_global_wgs, pt_global_wgs = parse_pbf_inside_aoi(str(pbf_path), aoi_poly)

    # Step 3 — Project to meters & prep Duomo
    banner("Step 3 — Projecting & prepping sources")
    utm_epsg = local_utm_epsg(center_lon, center_lat)
    crs_m = CRS.from_epsg(utm_epsg)

    to_m = Transformer.from_crs(4326, crs_m, always_xy=True).transform
    to_wgs = Transformer.from_crs(crs_m, 4326, always_xy=True).transform

    pois_m            = pois_wgs.to_crs(crs_m)            if len(pois_wgs)            else pois_wgs
    poc_m             = poc_wgs.to_crs(crs_m)             if len(poc_wgs)             else poc_wgs
    supers_global_m   = supers_global_wgs.to_crs(crs_m)   if len(supers_global_wgs)   else supers_global_wgs
    pt_global_m       = pt_global_wgs.to_crs(crs_m)       if len(pt_global_wgs)       else pt_global_wgs
    duomo_m_point     = shp_transform(to_m, Point(DUOMO_LON, DUOMO_LAT))

    # Build AEQD → UTM square that covers the AOI circle (for quadtree work in meters)
    aeqd = make_aeqd_crs(center_lon, center_lat)
    to_aeqd = Transformer.from_crs(4326, aeqd, always_xy=True).transform
    aoi_m   = shp_transform(to_aeqd, aoi_poly)
    to_utm  = Transformer.from_crs(aeqd, crs_m, always_xy=True).transform
    aoi_m_utm = shp_transform(to_utm, aoi_m)
    def bounding_square(poly_m):
        minx, miny, maxx, maxy = poly_m.bounds
        side = max(maxx-minx, maxy-miny)
        cx, cy = (minx+maxx)/2.0, (miny+maxy)/2.0
        return box(cx-side/2.0, cy-side/2.0, cx+side/2.0, cy+side/2.0)
    root_square_m = bounding_square(aoi_m_utm)

    # Step 4 — Quadtree
    leaves = build_quadtree(root_square_m, pois_m, min_tile_m=args.min_tile_m)

    # Step 5 — Distances (outside fallback), buffers, counts
    banner("Step 5 — Distances, buffers, and counts")
    records, poly_geoms_wgs = [], []
    iterator = enumerate(leaves, start=1)
    if TQDM: iterator = tqdm(iterator, total=len(leaves), desc="Tiles", unit="tile")

    # Spatial indexes (GLOBAL)
    pt_sidx    = pt_global_m.sindex     if len(pt_global_m)     else None
    super_sidx = supers_global_m.sindex if len(supers_global_m) else None

    # Spatial indexes (AOI) for counts inside tile
    poi_sidx = pois_m.sindex if len(pois_m) else None
    poc_sidx = poc_m.sindex  if len(poc_m)  else None

    def points_in_poly(poly_m: Polygon, points_gdf: gpd.GeoDataFrame, sidx):
        if sidx is None or len(points_gdf) == 0:
            return points_gdf.iloc[[]]
        cand = list(sidx.query(poly_m, predicate="intersects"))
        if not cand: return points_gdf.iloc[[]]
        sub = points_gdf.iloc[cand]
        mask = sub.within(gpd.GeoSeries([poly_m], crs=points_gdf.crs).iloc[0])
        return sub[mask]

    for i, (poly_m, depth, _) in iterator:
        # counts inside tile
        n_pois = int(len(points_in_poly(poly_m, pois_m, poi_sidx)))
        n_poc  = int(len(points_in_poly(poly_m, poc_m,  poc_sidx)))

        # Tile-level polygon distances (can be 0 if facility inside tile)
        if pt_sidx:
            d_pt_any = min_poly_to_nearest_point_any(poly_m, pt_global_m, pt_sidx)
            d_pt_out = min_poly_to_nearest_point_outside(poly_m, pt_global_m, pt_sidx)
        else:
            d_pt_any = float("nan"); d_pt_out = float("nan")

        if super_sidx:
            d_super_any = min_poly_to_nearest_point_any(poly_m, supers_global_m, super_sidx)
            d_super_out = min_poly_to_nearest_point_outside(poly_m, supers_global_m, super_sidx)
        else:
            d_super_any = float("nan"); d_super_out = float("nan")

        d_duomo_poly = float(poly_m.distance(duomo_m_point))

        # Fallbacks
        d_pt_final    = d_pt_out    if (math.isfinite(d_pt_any)    and d_pt_any    == 0.0 and math.isfinite(d_pt_out))    else d_pt_any
        d_super_final = d_super_out if (math.isfinite(d_super_any) and d_super_any == 0.0 and math.isfinite(d_super_out)) else d_super_any

        if math.isfinite(d_duomo_poly) and d_duomo_poly == 0.0:
            # Duomo inside tile → use tile-centroid distance
            minx, miny, maxx, maxy = poly_m.bounds
            cx, cy = (minx+maxx)/2.0, (miny+maxy)/2.0
            d_duomo_final = float(Point(cx, cy).distance(duomo_m_point))
        else:
            d_duomo_final = d_duomo_poly

        # buffer counts (compact names)
        buf_counts = {}
        for r in PT_BUFFER_LIST:
            buf_counts[f"pt{r}_cnt"] = count_points_within_buffer(poly_m, pt_global_m, pt_sidx, r)
        for r in SUPER_BUFFER_LIST:
            buf_counts[f"super{r}_cnt"] = count_points_within_buffer(poly_m, supers_global_m, super_sidx, r)

        # geometry for output (to WGS)
        minx, miny, maxx, maxy = poly_m.bounds
        cx, cy = (minx+maxx)/2.0, (miny+maxy)/2.0
        center_wgs = shp_transform(Transformer.from_crs(crs_m, 4326, always_xy=True).transform, Point(cx, cy))
        poly_wgs   = shp_transform(Transformer.from_crs(crs_m, 4326, always_xy=True).transform, poly_m)
        poly_geoms_wgs.append(poly_wgs)

        rec = {
            "tile_id": i,
            "qt_depth": depth,
            "n_pois": n_pois,
            "has_poi": 1 if n_pois > 0 else 0,
            "n_poc": n_poc,
            "has_poc": 1 if n_poc > 0 else 0,
            "lon": center_wgs.x,
            "lat": center_wgs.y,
            "side_m": maxx - minx,
            "area_m2": (maxx - minx) * (maxy - miny),

            # Lean tile-level distances with outside fallback
            "d_pt_m": d_pt_final,
            "d_super_m": d_super_final,
            "d_duomo_m": d_duomo_final,
        }
        rec.update(buf_counts)
        records.append(rec)

    # Step 6 — Write files (CSV + GPKG)
    banner("Step 6 — Writing files")
    csv_grid  = os.path.join(out_dir, "quadtree_grid.csv")   # same name, lean columns
    gpkg_path = os.path.join(out_dir, "quadtree_layers.gpkg")

    # CSV with a short commented preamble explaining distance semantics
    df_out = pd.DataFrame.from_records(records)
    preamble = [
        "# Column guide:\n",
        "# d_super_m : tile→nearest supermarket (m); if 0, uses nearest outside-of-tile supermarket.\n",
        "# d_pt_m    : tile→nearest public-transport stop (m); if 0, uses nearest outside-of-tile stop.\n",
        "# d_duomo_m : tile→Duomo (m); if Duomo inside the tile, uses tile-centroid→Duomo.\n",
        "# pt200_cnt/pt400_cnt, super300_cnt/super600_cnt : counts of facilities within the given radius from the tile polygon.\n",
        "\n"
    ]
    with open(csv_grid, "w", encoding="utf-8", newline="") as f:
        f.writelines(preamble)
        df_out.to_csv(f, index=False, encoding="utf-8")
    success(f"CSV written: {csv_grid}")

    # GPKG (for GIS visualization) — unchanged layers
    grid_gdf = gpd.GeoDataFrame(records, geometry=poly_geoms_wgs, crs="EPSG:4326")
    try:
        aoi_gdf = gpd.GeoDataFrame(geometry=[aoi_poly], crs="EPSG:4326")
        grid_gdf = gpd.overlay(grid_gdf, aoi_gdf, how="intersection", keep_geom_type=True)
    except Exception as e:
        warn(f"Overlay clip failed ({e}); writing unclipped grid.")

    # Ensure fresh GPKG
    if os.path.exists(gpkg_path):
        try: os.remove(gpkg_path)
        except Exception: pass

    grid_gdf.to_file(gpkg_path, layer="quadtree_grid", driver="GPKG")
    (fiona_safe(pois_wgs) if len(pois_wgs) else gpd.GeoDataFrame(pois_wgs, geometry="geometry", crs="EPSG:4326")) \
        .to_file(gpkg_path, layer="pois_target", driver="GPKG")
    (fiona_safe(poc_wgs) if len(poc_wgs) else gpd.GeoDataFrame(poc_wgs, geometry="geometry", crs="EPSG:4326")) \
        .to_file(gpkg_path, layer="pocs_control", driver="GPKG")
    (fiona_safe(supers_wgs) if len(supers_wgs) else gpd.GeoDataFrame(supers_wgs, geometry="geometry", crs="EPSG:4326")) \
        .to_file(gpkg_path, layer="supermarkets", driver="GPKG")
    (fiona_safe(pt_wgs) if len(pt_wgs) else gpd.GeoDataFrame(pt_wgs, geometry="geometry", crs="EPSG:4326")) \
        .to_file(gpkg_path, layer="transport", driver="GPKG")
    (fiona_safe(supers_global_wgs) if len(supers_global_wgs) else gpd.GeoDataFrame(supers_global_wgs, geometry="geometry", crs="EPSG:4326")) \
        .to_file(gpkg_path, layer="supermarkets_all", driver="GPKG")
    (fiona_safe(pt_global_wgs) if len(pt_global_wgs) else gpd.GeoDataFrame(pt_global_wgs, geometry="geometry", crs="EPSG:4326")) \
        .to_file(gpkg_path, layer="transport_all", driver="GPKG")
    success(f"GPKG written: {gpkg_path}")
    banner("All done")

if __name__ == "__main__":
    try:
        main()
    finally:
        # Ensure clean return to the calling shell (no lingering state)
        try:
            sys.stdout.flush(); sys.stderr.flush()
        except Exception:
            pass
        os._exit(0)
