@echo off
cd /d "%~dp0.."
python run_docking.py %*
echo.
pause
