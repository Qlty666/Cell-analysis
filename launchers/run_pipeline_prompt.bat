@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0.."
call "%~dp0_common.bat"
if errorlevel 1 (
  pause
  exit /b 1
)
set /p ACC=Input GSE accession e.g. GSE125449:
set /p OUT=Input output folder e.g. results\liver_cancer:
set /p SP=Input species: hs=human, mm=mouse, default hs:
if "%ACC%"=="" goto err_acc
if "%OUT%"=="" goto err_out
if "%SP%"=="" set SP=hs
%PYTHON_EXE% scripts\run_pipeline.py "%ACC%" --output "%OUT%" --species %SP%
set "RC=%errorlevel%"
echo.
pause
exit /b %RC%

:err_acc
echo ERROR: GSE accession is required.
pause
exit /b 1

:err_out
echo ERROR: output folder is required.
pause
exit /b 1
