@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
call "%~dp0launchers\_common.bat"
if errorlevel 1 (
  pause
  exit /b 1
)
%PYTHON_EXE% scripts\liverbio.py %*
exit /b %errorlevel%
