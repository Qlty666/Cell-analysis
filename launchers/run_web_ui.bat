@echo off
cd /d "%~dp0.."
echo Starting unified web UI. Default: /full full pipeline. Pages: / single-cell, /dock virtual screening, /results results manifest, /tasks task progress.
echo Use --page dock, --page full, --page tasks or --page results to open a specific page.
python web\web_ui.py %*
pause
