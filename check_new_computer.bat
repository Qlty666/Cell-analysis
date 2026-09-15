@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
call "%~dp0launchers\_common.bat"
if errorlevel 1 (
  pause
  exit /b 1
)
echo Checking the full environment...
%PYTHON_EXE% launchers\install_environment.py check full
set "RC=%errorlevel%"
if not "%RC%"=="0" (
  echo Environment check found missing components.
  pause
  exit /b %RC%
)
echo.
echo Environment check passed.
pause
exit /b %RC%
