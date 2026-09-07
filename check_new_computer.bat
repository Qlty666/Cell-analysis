@echo off
chcp 65001 >nul
cd /d "%~dp0"
echo Checking the full environment...
python launchers\install_environment.py check full
if errorlevel 1 (
  echo Environment check found missing components.
  pause
  exit /b 1
)
echo.
echo Environment check passed.
pause
exit /b 0
