@echo off
REM Convenience wrapper: latticegen2.bat -i part.step -cc 10 -t 1.5 [options]
REM
REM Picks an interpreter in this order:
REM   1. %LATTICEGEN2_PYTHON%          explicit override always wins
REM   2. .\runtime\python.exe          portable release bundle, nothing installed
REM   3. .\.venv\Scripts\python.exe    wheels release bundle, after install.bat
REM   4. python                        a plain checkout on a prepared machine
REM
REM The tool runs straight from this directory with no install step, so an
REM offline workstation needs either a release bundle or Python plus the two
REM dependencies listed in README.md.
setlocal
set "PY="
if not "%LATTICEGEN2_PYTHON%"=="" set "PY=%LATTICEGEN2_PYTHON%"
if "%PY%"=="" if exist "%~dp0runtime\python.exe" set "PY=%~dp0runtime\python.exe"
if "%PY%"=="" if exist "%~dp0.venv\Scripts\python.exe" set "PY=%~dp0.venv\Scripts\python.exe"
if "%PY%"=="" set "PY=python"
REM Always launch src/main.py rather than the installed console script: boundary
REM workers use multiprocessing "spawn", which re-imports this __main__ module in
REM each child, and main.py is what puts src/ on sys.path for them.
"%PY%" "%~dp0src\main.py" %*
exit /b %ERRORLEVEL%
