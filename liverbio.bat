@echo off
setlocal
cd /d "%~dp0"
python scripts\liverbio.py %*
exit /b %errorlevel%
