@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0.."
call "%~dp0_common.bat"
if errorlevel 1 (
  pause
  exit /b 1
)
%PYTHON_EXE% scripts\run_molecular_docking.py %*
set "RC=%errorlevel%"
echo.
pause
exit /b %RC%
