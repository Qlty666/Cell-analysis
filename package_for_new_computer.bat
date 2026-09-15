@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
call "%~dp0launchers\_common.bat"
if errorlevel 1 (
  pause
  exit /b 1
)
%PYTHON_EXE% launchers\package_portable.py %*
set "RC=%errorlevel%"
if not "%RC%"=="0" (
  echo Packaging failed.
  pause
  exit /b %RC%
)
pause
exit /b %RC%
