@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
call "%~dp0launchers\_common.bat"
if errorlevel 1 (
  pause
  exit /b 1
)
echo Installing the full environment for a new computer...
%PYTHON_EXE% launchers\install_environment.py install full --with-ml
set "RC=%errorlevel%"
if not "%RC%"=="0" (
  echo Setup failed. See the messages above.
  pause
  exit /b %RC%
)
echo.
echo Setup complete. Run check_new_computer.bat to verify, then use liverbio.bat.
pause
exit /b %RC%
