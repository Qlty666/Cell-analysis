@echo off
cd /d "%~dp0.."
python launchers\install_environment.py check expression %*
pause
