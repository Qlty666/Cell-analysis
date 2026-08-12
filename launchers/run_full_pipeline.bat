@echo off
setlocal
cd /d "%~dp0.."
python scripts\run_full_pipeline.py %*
exit /b %errorlevel%
