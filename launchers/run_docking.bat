@echo off
cd /d "%~dp0.."
python scripts\run_docking.py %*
echo.
pause
