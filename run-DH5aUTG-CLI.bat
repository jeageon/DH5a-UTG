@echo off
setlocal EnableExtensions EnableDelayedExpansion

set "PROJECT_DIR=%~dp0"

if exist "%PROJECT_DIR%.venv\Scripts\python.exe" (
  set "PY_CMD=%PROJECT_DIR%.venv\Scripts\python.exe"
  set "PY_ARGS="
) else (
  where py >nul 2>nul
  if errorlevel 1 (
    where python >nul 2>nul
    if errorlevel 1 (
      echo Python not found. Please install Python 3 and ensure PATH is configured.
      pause
      exit /b 1
    )
    set "PY_CMD=python"
    set "PY_ARGS="
  ) else (
    set "PY_CMD=py"
    set "PY_ARGS=-3"
  )
)

"%PY_CMD%" %PY_ARGS% -m pip install --upgrade pip >nul
"%PY_CMD%" %PY_ARGS% -m pip install -r "%PROJECT_DIR%requirements.txt"

cd /d "%PROJECT_DIR%"
"%PY_CMD%" %PY_ARGS% -m src.main %*
exit /b %ERRORLEVEL%

