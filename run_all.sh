#!/usr/bin/env bash
set -euo pipefail

# -------- Config (edit here or override via environment) --------
# Geofabrik: Nord-Ovest extract (Milan region)
PBF_URL="${PBF_URL:-https://download.geofabrik.de/europe/italy/nord-ovest-latest.osm.pbf}"
PBF_PATH="${PBF_PATH:-data/pbf/nord-ovest-latest.osm.pbf}"

# Area of interest (Duomo-ish center; meters radius)
CENTER="${CENTER:-45.4642,9.1914}"   # lat,lon
RADIUS_M="${RADIUS_M:-1500}"

# Output base; Python creates data/quadtree_YYYYMMDD_HHMMSS/
OUT_DIR="${OUT_DIR:-data}"
# ---------------------------------------------------------------

mkdir -p "$(dirname "$PBF_PATH")" "$OUT_DIR"

# --- Create/use local virtualenv & install Python deps ---
if [ ! -d ".venv" ]; then
  python3 -m venv .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt

# --- Auto-download the PBF if missing (kept under data/pbf/) ---
if [ ! -f "$PBF_PATH" ]; then
  echo "PBF missing — downloading:"
  echo "  $PBF_URL"
  python - <<'PY'
import os, ssl, urllib.request
url = os.environ["PBF_URL"]; dst = os.environ["PBF_PATH"]
os.makedirs(os.path.dirname(dst), exist_ok=True)
ctx = ssl.create_default_context()
with urllib.request.urlopen(url, context=ctx) as r, open(dst, "wb") as f:
    while True:
        b = r.read(1<<15)
        if not b: break
        f.write(b)
print("Download complete:", dst)
PY
fi

# --- Run the sampler (writes a timestamped run folder under data/) ---
python python/poi_point2point_distances.py \
  --center "$CENTER" --radius-m "$RADIUS_M" \
  --pbf "$PBF_PATH" --out-dir "$OUT_DIR" --min-tile-m 50

# --- Hint the newest CSV to pick in R (we do NOT pass it automatically) ---
CSV_HINT="$(ls -1t "$OUT_DIR"/quadtree_*/quadtree_grid.csv 2>/dev/null | head -n1 || true)"
if [ -n "$CSV_HINT" ]; then
  echo
  echo "When the R picker opens, choose this file:"
  echo "  $CSV_HINT"
  echo
fi

# --- Run the models (picker-only by design) ---
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Error: Rscript not found in PATH. Please install R (>= 4.2) and rerun." >&2
  exit 1
fi
Rscript R/run_all_proximity_models_single_txt.R

# --- Point to the final one-stop report ---
LATEST_REPORT="$(ls -1dt "$OUT_DIR"/quadtree_*/*/ 2>/dev/null | grep model_summaries | head -n1 || true)"
if [ -n "$LATEST_REPORT" ]; then
  echo "Done. See: ${LATEST_REPORT}proximity_models_all.txt"
else
  echo "Done. Look under the newest model_summaries_*/proximity_models_all.txt"
fi
