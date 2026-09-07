@echo off
cd /d "%~dp0.."
python launchers\install_environment.py check molecular-docking %*
pause
