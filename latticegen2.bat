@echo off
REM Convenience wrapper: latticegen2.bat -i part.step -cc 10 -t 1.5 [options]
REM
REM Uses the LATTICEGEN2_PYTHON environment variable if set, otherwise whatever
REM "python" resolves to. The tool runs straight from this checkout, with no
REM install step, so an offline workstation only needs Python plus the two
REM dependencies listed in README.md.
setlocal
if "%LATTICEGEN2_PYTHON%"=="" (set PY=python) else (set PY=%LATTICEGEN2_PYTHON%)
"%PY%" "%~dp0src\main.py" %*
exit /b %ERRORLEVEL%
