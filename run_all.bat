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

mkdir "%OUT_DIR%" >nul 2>nul
mkdir "%~dp0data\pbf" >nul 2>nul

REM --- Download PBF if missing (via curl) ---
if not exist "%PBF_PATH%" (
  echo PBF missing — downloading:
  echo   %PBF_URL%

  where curl >nul 2>nul
  if errorlevel 1 (
    echo Error: 'curl' command not found.
    echo Please install curl or download the file manually to:
    echo   %PBF_PATH%
    pause
    exit /b 1
  )
  
  curl -L -o "%PBF_PATH%" "%PBF_URL%"
  if errorlevel 1 (
      echo.
      echo Error: curl failed to download the file.
      echo Please check your internet connection and the URL.
      pause
      exit /b 1
  )
  echo Download complete.
)

REM --- Python venv & deps ---
echo Setting up Python environment...
if not exist .venv (
  python -m venv .venv
)
call .venv\Scripts\activate.bat
python -m pip install --upgrade pip >nul 
pip install -r requirements.txt

REM --- Run the sampler (writes a timestamped run folder under data\) ---
python "%PY_SCRIPT%" ^
  --center "%CENTER%" --radius-m %RADIUS_M% ^
  --pbf "%PBF_PATH%" --out-dir "%OUT_DIR%" --min-tile-m 50

REM --- Find the newest CSV to pass to R ---
for /f "usebackq delims=" %%i in (`powershell -NoProfile -Command ^
  "(Get-ChildItem -Path '%OUT_DIR%\quadtree_*\quadtree_grid.csv' -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1).FullName"`) do set "CSV_HINT=%%i"

if "%CSV_HINT%"=="" (
  echo Error: No quadtree_grid.csv file found for the R script to process.
  pause
  exit /b 1
)

REM --- Run the models ---
where Rscript >nul 2>nul
if errorlevel 1 (
  echo Error: Rscript not found in PATH. Please install R ^(^>= 4.2^) and rerun.
  pause
  exit /b 1
)
Rscript R\fullanalysis.r "%CSV_HINT%" 

REM --- Point to the final one-stop report ---
for /f "usebackq delims=" %%i in (`powershell -NoProfile -Command ^
  "(Get-ChildItem -Path '%OUT_DIR%\quadtree_*' -Directory -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | ForEach-Object { Get-ChildItem -Path $_.FullName -Directory | Where-Object { $_.Name -like 'model_summaries_*' } } | Select-Object -First 1).FullName"`) do set "LATEST=%%i"
if not "%LATEST%"=="" (
  echo Done. See: %LATEST%\proximity_models_all.txt
) else (
  echo Done. Look under the newest model_summaries_*\proximity_models_all.txt
)
pause
