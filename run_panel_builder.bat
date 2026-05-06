@echo off
REM Windows batch file for Panel Builder
REM Generates panel dataset from sampled establishments

setlocal enabledelayedexpansion

set "CONFIG_FILE=config.yml"
set "PYTHON_CMD=python"

echo.
echo ================================================================================
echo.
echo                    OSM PANEL BUILDER (Windows)
echo.
echo ================================================================================
echo.

REM Get script directory
cd /d "%~dp0"

REM Check for config file
if not exist "%CONFIG_FILE%" (
    echo. [FAIL] Config file not found: %CONFIG_FILE%
    pause
    exit /b 1
)

echo. [OK] Using config: %CONFIG_FILE%

REM Get the sample directory (most recent)
if "%~1"=="" (
    echo.
    echo. [>>] Config will auto-detect most recent sample
) else (
    set "SAMPLE_DIR=%~1"
    if not exist "!SAMPLE_DIR!" (
        echo. [FAIL] Directory not found: !SAMPLE_DIR!
        pause
        exit /b 1
    )
    set "INPUT_FILE=!SAMPLE_DIR!\establishments.csv"
    if not exist "!INPUT_FILE!" (
        echo. [FAIL] Establishments file not found: !INPUT_FILE!
        pause
        exit /b 1
    )
    echo.
    echo. [>>] Using specified sample: !INPUT_FILE!
    echo. [INPUT_ARG]=--input "!INPUT_FILE!"
)

echo.
echo --------------------------------------------------------------------------------
echo. STEP 1: Checking Dependencies
echo --------------------------------------------------------------------------------
echo.

echo Checking Python...
%PYTHON_CMD% --version >nul 2>&1
if errorlevel 1 (
    echo. [FAIL] Python not found
    pause
    exit /b 1
)
echo. [OK] Python available

echo Checking Python packages...
%PYTHON_CMD% -c "import geopandas; import pandas; import yaml" >nul 2>&1
if errorlevel 1 (
    echo. [FAIL] Required Python packages not found
    echo. Install with: pip install geopandas pandas pyyaml
    pause
    exit /b 1
)
echo. [OK] Packages available

echo.
echo --------------------------------------------------------------------------------
echo. STEP 2: Generate Panel
echo --------------------------------------------------------------------------------
echo.

REM Run panel builder
if defined INPUT_ARG (
    %PYTHON_CMD% scripts\panel_builder.py --config "%CONFIG_FILE%" %INPUT_ARG%
) else (
    %PYTHON_CMD% scripts\panel_builder.py --config "%CONFIG_FILE%"
)

if errorlevel 1 (
    echo.
    echo ================================================================================
    echo.
    echo. [FAIL] Panel generation FAILED
    echo.
    echo ================================================================================
    echo.
    echo. Check the error messages above
    echo.
    pause
    exit /b 1
)

echo.
echo ================================================================================
echo.
echo. [OK] SUCCESS - Panel dataset generated successfully!
echo.
echo ================================================================================
echo.
echo. Check data\output\panels\ for the results
echo.
pause
exit /b 0
