@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0.."
call "%~dp0_common.bat"
if errorlevel 1 (
  pause
  exit /b 1
)
echo Select the module environment to install:
echo   1. Expression analysis
echo   2. Dataset search
echo   3. Virtual screening / docking
echo   4. Standalone molecular docking
echo   5. Molecular dynamics (GROMACS)
echo   6. Full integrated pipeline
echo   7. Web console
echo   8. Project Codex skills
set /p env_choice=Enter number:
set "MOD="
if "%env_choice%"=="1" set "MOD=expression"
if "%env_choice%"=="2" set "MOD=datasets"
if "%env_choice%"=="3" set "MOD=docking"
if "%env_choice%"=="4" set "MOD=molecular-docking"
if "%env_choice%"=="5" set "MOD=md"
if "%env_choice%"=="6" set "MOD=full"
if "%env_choice%"=="7" set "MOD=web"
if "%env_choice%"=="8" set "MOD=skills"
if not defined MOD (
  echo ERROR: invalid selection "%env_choice%".
  pause
  exit /b 2
)
%PYTHON_EXE% launchers\install_environment.py install %MOD%
set "RC=%errorlevel%"
pause
exit /b %RC%
