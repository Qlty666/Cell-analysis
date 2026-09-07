@echo off
chcp 65001 >nul
cd /d "%~dp0"
echo Installing the full environment for a new computer...
python launchers\install_environment.py install full --with-ml
if errorlevel 1 (
  echo Setup failed. See the messages above.
  pause
  exit /b 1
)
echo.
echo Setup complete. Run check_new_computer.bat to verify, then use liverbio.bat.
pause
exit /b 0
