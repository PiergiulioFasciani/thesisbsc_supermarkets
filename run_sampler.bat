@echo off
REM Windows batch file for OSM Sampler
REM Creates/activates conda environment and runs the sampler

setlocal enabledelayedexpansion

REM Color codes and symbols (basic Windows console)
set "CHECK=[OK]"
set "CROSS=[FAIL]"
set "ARROW=[>>]"
set "INFO=[*]"
set "WARN=[!]"

echo.
echo ================================================================================
echo.
echo               OSM ESTABLISHMENT SAMPLER - SETUP ^& RUN (Windows)
echo.
echo ================================================================================
echo.

REM Get script directory
cd /d "%~dp0"

REM Configuration
set "ENV_NAME=osm_sampler_env"
set "CONFIG_FILE=config.yml"
set "PYTHON_CMD=python"

echo [STEP 1/4] Checking Dependencies
echo ================================================================================
echo.

REM Check if Python is available
echo Checking for Python...
%PYTHON_CMD% --version >nul 2>&1
if errorlevel 1 (
    echo %CROSS%
    echo ERROR: Python is not installed or not in PATH
    echo Install from: https://www.python.org/downloads/
    pause
    exit /b 1
)
echo %CHECK%

echo.
echo [STEP 2/4] Validating Configuration
echo ================================================================================
echo.

REM Check if config file exists
echo Checking config file (%CONFIG_FILE%)...
if not exist "%CONFIG_FILE%" (
    echo %CROSS%
    echo ERROR: Config file not found
    pause
    exit /b 1
)
echo %CHECK%

REM Check if input file exists
for /f "delims=" %%A in (
    '%PYTHON_CMD% -c "import yaml; cfg=yaml.safe_load(open('%CONFIG_FILE%')); s=cfg.get('sampler', cfg); print(s.get('input_file', ''))"'
) do set "INPUT_FILE=%%A"

echo Checking input file...
if not exist "%INPUT_FILE%" (
    echo %WARN%
    echo.
    echo WARNING: Input file not found: %INPUT_FILE%
    echo Place your OSM history file (.osh.pbf) in: data\input\
    echo Or update 'input_file' in the sampler section of %CONFIG_FILE%
    echo.
    setlocal disabledelayedexpansion
    set /p "CONTINUE=Continue anyway? [y/N] "
    setlocal enabledelayedexpansion
    if /i not "%CONTINUE%"=="y" (
        echo %CROSS% Aborted
        exit /b 1
    )
) else (
    echo %CHECK%
    for %%A in ("%INPUT_FILE%") do set "FILE_SIZE=%%~zA bytes"
    echo %ARROW% File size: !FILE_SIZE!
)

REM Display config info
echo.
echo Configuration Summary:
for /f "delims=" %%A in (
    '%PYTHON_CMD% -c "import yaml; cfg=yaml.safe_load(open('%CONFIG_FILE%')); s=cfg.get('sampler', cfg); print(s.get('center_lat', ''))"'
) do set "CENTER_LAT=%%A"

for /f "delims=" %%A in (
    '%PYTHON_CMD% -c "import yaml; cfg=yaml.safe_load(open('%CONFIG_FILE%')); s=cfg.get('sampler', cfg); print(s.get('center_lon', ''))"'
) do set "CENTER_LON=%%A"

for /f "delims=" %%A in (
    '%PYTHON_CMD% -c "import yaml; cfg=yaml.safe_load(open('%CONFIG_FILE%')); s=cfg.get('sampler', cfg); print(s.get('radius_meters', ''))"'
) do set "RADIUS=%%A"

for /f "delims=" %%A in (
    '%PYTHON_CMD% -c "import yaml; cfg=yaml.safe_load(open('%CONFIG_FILE%')); s=cfg.get('sampler', cfg); print(s.get('start_date', ''))"'
) do set "START_DATE=%%A"

for /f "delims=" %%A in (
    '%PYTHON_CMD% -c "import yaml; cfg=yaml.safe_load(open('%CONFIG_FILE%')); s=cfg.get('sampler', cfg); print(s.get('end_date', ''))"'
) do set "END_DATE=%%A"

echo %ARROW% Center: %CENTER_LAT%, %CENTER_LON%
echo %ARROW% Radius: %RADIUS%m
echo %ARROW% Date range: %START_DATE% to %END_DATE%

echo.
echo [STEP 3/4] Verifying Dependencies
echo ================================================================================
echo.

echo Checking Python packages...
%PYTHON_CMD% -c "import osmium; import geopandas; import yaml; import shapely; import pyproj" >nul 2>&1
if errorlevel 1 (
    echo %CROSS%
    echo ERROR: Package verification failed
    echo Please ensure the required packages are installed:
    echo   pip install osmium geopandas pyyaml shapely pyproj tqdm
    pause
    exit /b 1
)
echo %CHECK%

echo.
echo [STEP 4/4] Running Sampler
echo ================================================================================
echo.

%PYTHON_CMD% scripts\osm_sampler.py --config "%CONFIG_FILE%"

if errorlevel 1 (
    echo.
    echo ================================================================================
    echo.
    echo %CROSS% SAMPLING FAILED
    echo.
    echo ================================================================================
    echo.
    echo Check the error messages above for details
    echo.
    pause
    exit /b 1
)

echo.
echo ================================================================================
echo.
echo %CHECK% SAMPLING COMPLETED SUCCESSFULLY
echo.
echo ================================================================================
echo.
echo Next steps:
echo %ARROW% Check output in: data\output\samples\
echo %ARROW% Open GeoPackage in QGIS for visualization
echo %ARROW% Analyze CSV with pandas/R
echo.
echo To run again:
echo %ARROW% run_sampler.bat
echo.
pause
exit /b 0
