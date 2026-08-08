@echo off
cd /d "%~dp0.."
echo Starting unified web UI. Pages: / single-cell, /dock virtual screening, /full full pipeline.
echo Use --page dock or --page full to open a specific page.
python web\web_ui.py %*
pause
