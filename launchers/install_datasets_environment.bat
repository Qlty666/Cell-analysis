@echo off
cd /d "%~dp0.."
python launchers\install_environment.py install datasets %*
pause
