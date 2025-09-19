#!/usr/bin/env bash
set -euo pipefail

# -------- Config (edit here or override via env) --------
# Geofabrik Nord-Ovest (Milan region)
PBF_URL="${PBF_URL:-https://download.geofabrik.de/europe/italy/nord-ovest-latest.osm.pbf}"
PBF_PATH="${PBF_PATH:-data/pbf/nord-ovest-latest.osm.pbf}"

# Area of interest (Duomo-ish)
CENTER="${CENTER:-45.4642,9.1914}"   # lat,lon
RADIUS_M="${RADIUS_M:-1500}"         # meters

# Output base; Python creates data/quadtree_YYYYMMDD_HHMMSS/
OUT_DIR="${OUT_DIR:-data}"

# Which Python entry? default main.py; set USE_LEGACY=1 to run legacy_main.py
PY_SCRIPT="python/main.py"
if [ "${USE_LEGACY:-0}" = "1" ]; then PY_SCRIPT="python/legacy_main.py"; fi
# -------------------------------------------------------

mkdir -p "$(dirname "$PBF_PATH")" "$OUT_DIR"

# --- Python venv & deps ---
if [ ! -d ".venv" ]; then python3 -m venv .venv; fi
# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt

# --- Download PBF if missing ---
if [ ! -f "$PBF_PATH" ]; then
  echo "PBF missing — downloading:"
  echo "  $PBF_URL"
  # Use curl with -L to follow redirects and --create-dirs to make the directory
  curl -L "$PBF_URL" -o "$PBF_PATH" --create-dirs
  echo "Download complete: $PBF_PATH"
fi

# --- Run the sampler (writes a timestamped run folder under data/) ---
python "$PY_SCRIPT" \
  --center "$CENTER" --radius-m "$RADIUS_M" \
  --pbf "$PBF_PATH" --out-dir "$OUT_DIR" --min-tile-m 50

# --- Find the newest CSV to pass to R ---
CSV_HINT="$(ls -1t "$OUT_DIR"/quadtree_*/quadtree_grid.csv 2>/dev/null | head -n1 || true)"
if [ -z "$CSV_HINT" ]; then
  echo "Error: No quadtree_grid.csv file found for the R script to process." >&2
  exit 1
fi

# --- Run the models ---
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Error: Rscript not found in PATH. Please install R (>= 4.2) and rerun." >&2
  exit 1
fi
echo "Running R analysis on: $CSV_HINT"
Rscript R/fullanalysis.r "$CSV_HINT"

# --- Point to the final one-stop report ---
LATEST_REPORT="$(ls -1dt "$OUT_DIR"/quadtree_*/*/ 2>/dev/null | grep model_summaries | head -n1 || true)"
if [ -n "$LATEST_REPORT" ]; then
  echo "Done. See: ${LATEST_REPORT}proximity_models_all.txt"
else
  echo "Done. Look under the newest model_summaries_*/proximity_models_all.txt"
fi
