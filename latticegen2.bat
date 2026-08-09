@echo off
REM Convenience wrapper for Windows. See specification.md §2 (Packaging form)
REM and README.md for usage. Forwards all arguments to src/main.jl using the
REM Julia project defined by Project.toml/Manifest.toml in this directory.
REM
REM `-t auto` enables Julia's own thread pool (distinct from the --cores/
REM --workers Distributed *process* pool used for tiling, docs/algorithm.md
REM §6.2): the classification stage's per-strut loop is threaded
REM (docs/algorithm.md §11.2) and is otherwise silently single-threaded, since
REM Julia defaults to 1 thread unless told otherwise.
setlocal
set SCRIPT_DIR=%~dp0
julia -t auto --project="%SCRIPT_DIR%" "%SCRIPT_DIR%src\main.jl" %*
exit /b %ERRORLEVEL%
