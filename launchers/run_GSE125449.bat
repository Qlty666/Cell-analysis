@echo off
rem Backward-compatible alias: README.md and docs/project_structure.md still
rem reference this file name. The interactive prompts live in
rem run_pipeline_prompt.bat.
call "%~dp0run_pipeline_prompt.bat" %*
exit /b %errorlevel%
