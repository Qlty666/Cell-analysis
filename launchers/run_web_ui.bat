@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0.."
call "%~dp0_common.bat"
if errorlevel 1 (
  pause
  exit /b 1
)
echo Starting unified web UI. Default: /full full pipeline. CADD pages: /dock virtual screening, /md-simulation molecular dynamics, /knockout virtual knockout, /network network toxicology, /faers FAERS, /validation real-data validation.
echo Other pages: / single-cell, /datasets dataset search, /molecular-docking molecular docking, /results results manifest, /tasks task progress, /guide usage tutorial, /environment environment board.
echo Use --page dock, --page md-simulation, --page knockout, --page network, --page faers, --page validation, --page molecular-docking, --page full, --page tasks, --page datasets, --page guide or --page environment to open a specific page.
%PYTHON_EXE% web\web_ui.py %*
set "RC=%errorlevel%"
if not "%RC%"=="0" (
  echo.
  echo Web UI exited with an error. Press any key to close this window.
  pause > nul
)
exit /b %RC%
