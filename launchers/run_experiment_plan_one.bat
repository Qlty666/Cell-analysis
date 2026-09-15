@echo off
setlocal
call "%~dp0_common.bat"
if errorlevel 1 exit /b %errorlevel%
python "%~dp0..\scripts\run_experiment_plan_one.py" ^
  --config "%~dp0..\config\experiment_plan_one.json" ^
  --output-root "%~dp0..\experiment_plan_one_results" %*
exit /b %errorlevel%
