@echo off
chcp 65001 >nul
cd /d "%~dp0"
python launchers\package_portable.py %*
if errorlevel 1 (
  echo Packaging failed.
  pause
  exit /b 1
)
pause
exit /b 0
