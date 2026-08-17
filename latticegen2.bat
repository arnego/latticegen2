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
set "RC=%ERRORLEVEL%"
if "%RC%"=="0" exit /b 0

REM A non-zero code is normally the tool reporting a real failure, so say nothing
REM and pass it through. But it is also what an interpreter that cannot load the
REM dependencies produces *before* main.py runs at all, and that case cannot be
REM reported from inside Python: a native library aborting during import kills
REM the process outright -- MKL's "cannot load mkl_intel_thread.2.dll" is not
REM even catchable as an exception -- and it exits 2, colliding with the
REM parameter-error code.
REM
REM So probe, but only for codes that could actually be that. 3, 4, 5, 6 and 130
REM come from the pipeline's own error paths (docs/algorithm.md S10), so by then
REM the tool has demonstrably run and a probe is pure delay -- which after a
REM Ctrl+C would mean an unexplained pause following a shutdown the spec asks to
REM be clean. Everything else can be an interpreter that never started: 1, 2, and
REM the shell's own 9009 for a command it could not find.
for %%C in (3 4 5 6 130) do if "%RC%"=="%%C" exit /b %RC%

REM `%ERRORLEVEL% NEQ 0` rather than `not errorlevel 1`: `if errorlevel N` tests
REM >= N, so a probe that dies from a native fault (0xC0000005 arrives as a
REM negative code) would read as success and suppress the diagnostic.
"%PY%" -c "import numpy, OCP.TopoDS" >nul 2>&1
if %ERRORLEVEL% EQU 0 exit /b %RC%

REM `>&2 echo text` rather than `echo text >&2`: the redirect has to lead, or the
REM space before it is echoed as trailing whitespace on every line.
>&2 echo.
>&2 echo FAILED: the interpreter below cannot load latticegen2's dependencies, so
>&2 echo         the tool never ran and the exit code above is not its own.
>&2 echo.
>&2 echo   interpreter: %PY%
>&2 echo.
"%PY%" -c "import numpy, OCP.TopoDS" 1>&2
>&2 echo.
>&2 echo   Point LATTICEGEN2_PYTHON at a Python 3.11+ that has numpy and
>&2 echo   cadquery-ocp:
>&2 echo     set LATTICEGEN2_PYTHON=C:\path\to\python.exe
>&2 echo.
>&2 echo   See README.md "Installation". A release bundle carries its own interpreter
>&2 echo   and needs none of this.
exit /b 1
