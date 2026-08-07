@echo off
cd /d "%~dp0.."
echo Starting unified web UI (single-cell + docking). Use --page dock to open docking page.
python web\web_ui.py %*
pause
