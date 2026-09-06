@echo off
cd /d "%~dp0.."
echo Starting unified web UI. Default: /full full pipeline. CADD pages: /dock virtual screening, /md-simulation molecular dynamics, /knockout virtual knockout, /network network toxicology, /faers FAERS, /validation real-data validation.
echo Other pages: / single-cell, /datasets dataset search, /molecular-docking molecular docking, /results results manifest, /tasks task progress.
echo Use --page dock, --page md-simulation, --page knockout, --page network, --page faers, --page validation, --page molecular-docking, --page full, --page tasks, --page datasets or --page results to open a specific page.
python web\web_ui.py %*
if errorlevel 1 (
  echo.
  echo Web UI exited with an error. Press any key to close this window.
  pause > nul
)
