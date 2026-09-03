@echo off
cd /d "%~dp0.."
echo Starting unified web UI. Default: /full full pipeline. Pages: / single-cell, /datasets dataset search, /dock virtual screening, /molecular-docking molecular docking, /results results manifest, /tasks task progress.
echo Use --page dock, --page molecular-docking, --page full, --page tasks, --page datasets or --page results to open a specific page.
python web\web_ui.py %*
if errorlevel 1 (
  echo.
  echo Web UI exited with an error. Press any key to close this window.
  pause > nul
)
