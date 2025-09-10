@echo off
setlocal enabledelayedexpansion

REM -------- Config (edit here or set as env vars) --------
REM Geofabrik Nord-Ovest (Milan region)
if "%PBF_URL%"==""   set "PBF_URL=https://download.geofabrik.de/europe/italy/nord-ovest-latest.osm.pbf"
if "%PBF_PATH%"==""  set "PBF_PATH=data\pbf\nord-ovest-latest.osm.pbf"

REM Area of interest (Duomo-ish)
if "%CENTER%"==""    set "CENTER=45.4642,9.1914"
if "%RADIUS_M%"==""  set "RADIUS_M=1500"

REM Output base; Python creates data\quadtree_YYYYMMDD_HHMMSS\
if "%OUT_DIR%"==""   set "OUT_DIR=data"

REM Which Python entry? default main.py; set USE_LEGACY=1 to run legacy_main.py
if "%USE_LEGACY%"=="1" (set "PY_SCRIPT=python\legacy_main.py") else (set "PY_SCRIPT=python\main.py")
REM -------------------------------------------------------

if not exist "%OUT_DIR%" mkdir "%OUT_DIR%"
if not exist "%~dp0data\pbf" mkdir "%~dp0data\pbf"

REM --- Download PBF if missing (via Python) ---
if not exist "%PBF_PATH%" (
  echo PBF missing — downloading:
  echo   %PBF_URL%
  py - <<PY
import os, ssl, urllib.request
url=os.environ['PBF_URL']; dst=os.environ['PBF_PATH']
os.makedirs(os.path.dirname(dst), exist_ok=True)
ctx=ssl.create_default_context()
with urllib.request.urlopen(url, context=ctx) as r, open(dst,'wb') as f:
    while True:
        b=r.read(1<<15)
        if not b: break
        f.write(b)
print("Download complete:", dst)
PY
)

REM --- Python venv & deps ---
if not exist .venv (
  py -3 -m venv .venv
)
call .venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt

REM --- Run the sampler (writes a timestamped run folder under data\) ---
python "%PY_SCRIPT%" ^
  --center "%CENTER%" --radius-m %RADIUS_M% ^
  --pbf "%PBF_PATH%" --out-dir "%OUT_DIR%" --min-tile-m 50

REM --- Hint the newest CSV to pick in R (we do NOT pass it automatically) ---
for /f "usebackq delims=" %%i in (`powershell -NoProfile -Command ^
  "(Get-ChildItem -Path '%OUT_DIR%\quadtree_*\\quadtree_grid.csv' -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1).FullName"`) do set "CSV_HINT=%%i"
if not "%CSV_HINT%"=="" (
  echo.
  echo When the R picker opens, choose this file:
  echo   %CSV_HINT%
  echo.
)

REM --- Run the models (picker-only by design) ---
where Rscript >nul 2>nul
if errorlevel 1 (
  echo Error: Rscript not found in PATH. Please install R (>= 4.2) and rerun.
  pause
  exit /b 1
)
Rscript R\fullanalysis.r

REM --- Point to the final one-stop report ---
for /f "usebackq delims=" %%i in (`powershell -NoProfile -Command ^
  "(Get-ChildItem -Path '%OUT_DIR%\quadtree_*' -Directory -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | ForEach-Object { Get-ChildItem -Path $_.FullName -Directory | Where-Object { $_.Name -like 'model_summaries_*' } } | Select-Object -First 1).FullName"`) do set "LATEST=%%i"
if not "%LATEST%"=="" (
  echo Done. See: %LATEST%\proximity_models_all.txt
) else (
  echo Done. Look under the newest model_summaries_*\proximity_models_all.txt
)
pause
