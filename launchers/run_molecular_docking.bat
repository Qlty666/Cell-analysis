@echo off
cd /d "%~dp0.."
python scripts\run_molecular_docking.py %*
echo.
pause
