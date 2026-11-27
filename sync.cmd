@echo off
setlocal
cd /d "%~dp0"

powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0sync.ps1" %*

if errorlevel 1 (
    echo.
    echo Sync failed. See the message above.
    pause
)

endlocal
