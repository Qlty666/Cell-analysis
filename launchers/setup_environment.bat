@echo off
chcp 65001 >nul
cd /d "%~dp0.."
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
if "%env_choice%"=="1" python launchers\install_environment.py install expression
if "%env_choice%"=="2" python launchers\install_environment.py install datasets
if "%env_choice%"=="3" python launchers\install_environment.py install docking
if "%env_choice%"=="4" python launchers\install_environment.py install molecular-docking
if "%env_choice%"=="5" python launchers\install_environment.py install md
if "%env_choice%"=="6" python launchers\install_environment.py install full
if "%env_choice%"=="7" python launchers\install_environment.py install web
if "%env_choice%"=="8" python launchers\install_environment.py install skills
pause
