@echo off
REM Windows batch file for Stacked PPML Event-Study Analysis
REM Checks dependencies and runs the stacked PPML analysis

setlocal enabledelayedexpansion

set "CONFIG_FILE=config.yml"
set "PANEL_DIR=data\output\panels"

echo.
echo ===================================================================
echo.
echo   STACKED PPML EVENT-STUDY ANALYSIS - DEPENDENCY CHECK (Windows)
echo.
echo ===================================================================
echo.
echo Date: %date% %time%
echo.

REM ==================================================================
REM 1. Check R Installation
REM ==================================================================

echo -------------------------------------------------------------------
echo 1. Checking R installation...
echo -------------------------------------------------------------------
echo.

where Rscript >nul 2>&1
if errorlevel 1 (
    echo [FAIL] R is not installed or not in PATH
    echo Please install R from: https://cran.r-project.org/
    pause
    exit /b 1
)

for /f "delims=" %%A in ('Rscript --version 2^>^&1 ^| findstr /r ".*"') do set "R_VERSION=%%A"
echo [OK] R found: %R_VERSION%
echo.

REM ==================================================================
REM 2. Check Required R Packages
REM ==================================================================

echo -------------------------------------------------------------------
echo 2. Checking required R packages...
echo -------------------------------------------------------------------
echo.

setlocal enabledelayedexpansion
set "PACKAGES=fixest ggplot2 dplyr tidyr data.table yaml zoo car HonestDiD"
set "MISSING_COUNT=0"

for %%P in (%PACKAGES%) do (
    Rscript -e "if (!require('%%P', quietly = TRUE)) quit(status = 1)" >nul 2>&1
    if errorlevel 1 (
        echo   [FAIL] %%P (missing)
        set /a "MISSING_COUNT=!MISSING_COUNT! + 1"
        set "MISSING_PACKAGES=!MISSING_PACKAGES! %%P"
    ) else (
        echo   [OK] %%P
    )
)

echo.

REM ==================================================================
REM 3. Install Missing Packages
REM ==================================================================

if !MISSING_COUNT! gtr 0 (
    echo -------------------------------------------------------------------
    echo 3. Installing !MISSING_COUNT! missing package(s)...
    echo -------------------------------------------------------------------
    echo.

    for %%P in (!MISSING_PACKAGES!) do (
        echo Installing %%P...
        
        if "%%P"=="HonestDiD" (
            echo   Installing HonestDiD from GitHub (requires devtools)...
            Rscript -e "if (!require('devtools', quietly = TRUE)) { install.packages('devtools', repos = 'https://cloud.r-project.org/') }; devtools::install_github('asheshrambachan/HonestDiD')"
        ) else (
            Rscript -e "install.packages('%%P', repos = 'https://cloud.r-project.org/')"
        )

        if errorlevel 1 (
            echo   [FAIL] ERROR: Failed to install %%P
            pause
            exit /b 1
        ) else (
            echo   [OK] %%P installed successfully
        )
    )
    echo.
    echo All missing packages installed successfully!
) else (
    echo All required packages are already installed!
)

echo.

REM ==================================================================
REM 4. Verify Configuration File
REM ==================================================================

echo -------------------------------------------------------------------
echo 4. Verifying configuration...
echo -------------------------------------------------------------------
echo.

if not exist "%CONFIG_FILE%" (
    echo [FAIL] Configuration file not found: %CONFIG_FILE%
    echo Please create the configuration file before running analysis
    pause
    exit /b 1
)

echo [OK] Configuration file found: %CONFIG_FILE%
echo.
echo Configuration settings:
echo -------------------------------------------------------------------
for /f "tokens=*" %%A in ('findstr /e "time_aggregation estimator start_date end_date pre_periods post_periods cluster_var" "%CONFIG_FILE%"') do (
    echo   %%A
)
echo -------------------------------------------------------------------
echo.

REM ==================================================================
REM 5. Check for Panel Data
REM ==================================================================

echo -------------------------------------------------------------------
echo 5. Checking for panel data...
echo -------------------------------------------------------------------
echo.

if not exist "%PANEL_DIR%" (
    echo [FAIL] Panel directory not found: %PANEL_DIR%
    echo Please run panel_builder.py first to create panel data
    pause
    exit /b 1
)

setlocal enabledelayedexpansion
set "PANEL_COUNT=0"
for /f "delims=" %%A in ('dir /s /b "%PANEL_DIR%\panel.csv" 2^>nul ^| find /c /v ""') do (
    set "PANEL_COUNT=%%A"
)

if !PANEL_COUNT! equ 0 (
    echo [FAIL] No panel.csv files found in %PANEL_DIR%
    echo Please run panel_builder.py first to create panel data
    pause
    exit /b 1
)

echo [OK] Found !PANEL_COUNT! panel file(s)

for /f "delims=" %%A in ('python -c "from pathlib import Path; panel_files = [p for p in Path('data/output/panels').glob('panel_*/panel.csv') if p.is_file()]; latest = max(panel_files, key=lambda p: p.stat().st_mtime) if panel_files else None; print(latest if latest else 'None')"') do (
    set "LATEST_PANEL=%%A"
)

echo   Latest panel: !LATEST_PANEL!
echo.

REM ==================================================================
REM 6. Run Stacked PPML Event-Study Analysis
REM ==================================================================

echo ===================================================================
echo.
echo   RUNNING STACKED PPML EVENT-STUDY ANALYSIS
echo.
echo ===================================================================
echo.
echo This will perform:
echo   - Standard stacked PPML event studies for count, entries, exits
echo   - Wald pretrend tests on the closest pre-treatment periods
echo   - Support-weighted post-treatment PPML ATT summaries
echo   - HonestDiD sensitivity analysis for all specifications
echo.
echo Starting analysis...
echo.

set "START_TIME=%date% %time%"

Rscript scripts\stacked_ppml_event_study.R

if errorlevel 1 (
    echo.
    echo ===================================================================
    echo.
    echo   [FAIL] ANALYSIS FAILED
    echo.
    echo ===================================================================
    echo.
    echo Please check the error messages above for details
    echo.
    pause
    exit /b 1
)

echo.
echo ===================================================================
echo.
echo   [OK] ANALYSIS COMPLETED SUCCESSFULLY
echo.
echo ===================================================================
echo.

REM Find latest output directory
for /f "delims=" %%A in ('dir /b /ad "data\output\analysis" ^| sort /r ^| findstr "stacked_ppml_es"') do (
    set "OUTPUT_DIR=data\output\analysis\%%A"
    goto :found_output
)

:found_output
if defined OUTPUT_DIR (
    echo Results location:
    echo   Directory: !OUTPUT_DIR!
    echo.
    echo Output files:
    echo   - stacked_ppml_event_study_results.txt (Full results text file^)
    echo   - csv\complete_summary.csv             (Summary table^)
    echo   - plots\standard\                     (Event-study and HonestDiD PNGs^)
    echo.
    
    for /f "delims=" %%A in ('dir /s /b "!OUTPUT_DIR!\plots\*.png" 2^>nul ^| find /c /v ""') do (
        echo   Generated %%A plots
    )
    
    echo.
    if exist "!OUTPUT_DIR!\csv\complete_summary.csv" (
        echo Quick summary (see full results in stacked_ppml_event_study_results.txt^):
        echo -------------------------------------------------------------------
        REM Show first 10 lines of summary
        setlocal enabledelayedexpansion
        set "LINE_COUNT=0"
        for /f "delims=" %%A in ('type "!OUTPUT_DIR!\csv\complete_summary.csv"') do (
            set /a "LINE_COUNT=!LINE_COUNT! + 1"
            if !LINE_COUNT! leq 10 (
                echo %%A
            )
        )
        echo -------------------------------------------------------------------
    )
)

echo.
echo ===================================================================
echo.
pause
exit /b 0
