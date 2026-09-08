@echo off
rem Shared Python 3 detection for the launcher scripts.
rem Sets PYTHON_EXE to "py -3" or "python"; returns 0 on success, 1 when no
rem working Python 3 interpreter can be found.
setlocal
set "PYTHON_EXE="
where py >nul 2>nul
if not errorlevel 1 (
  py -3 -c "import sys" >nul 2>nul
  if not errorlevel 1 set "PYTHON_EXE=py -3"
)
if not defined PYTHON_EXE (
  where python >nul 2>nul
  if not errorlevel 1 (
    python -c "import sys" >nul 2>nul
    if not errorlevel 1 set "PYTHON_EXE=python"
  )
)
if not defined PYTHON_EXE (
  echo ERROR: Python 3 was not found.
  echo Install Python 3.10 or newer from https://www.python.org/downloads/
  echo and tick "Add python.exe to PATH", then run this script again.
  endlocal
  exit /b 1
)
endlocal & set "PYTHON_EXE=%PYTHON_EXE%"
exit /b 0
