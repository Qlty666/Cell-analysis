@echo off
setlocal
cd /d "%~dp0.."

if "%~1"=="" (
  set /p "RUN=Enter dataset accession or run directory: "
  python scripts\liverbio.py analysis-export %RUN%
) else (
  python scripts\liverbio.py analysis-export %*
)
exit /b %errorlevel%
