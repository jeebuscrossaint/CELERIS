@echo off
REM ============================================================================
REM CELERIS one-shot verify: configure -> build (CLI + GUI) -> run physics
REM selftest. This is the "Verify" beat of the work loop. Run from repo root:
REM   scripts\check.bat
REM Exit code 0 = built and selftest ran; nonzero = build/configure failed.
REM (Eyeball the selftest numbers: energy=1.000000, 2D vs grcwa |d|<1e-2, Strehl ~0.55.)
REM ============================================================================
setlocal
cd /d "%~dp0.."
set "CMAKE=%USERPROFILE%\scoop\shims\cmake.exe"

echo [1/3] configure...
"%CMAKE%" -B build-msvc -DCELERIS_BUILD_GUI=ON >nul 2>&1 || (echo   CONFIGURE FAILED & exit /b 1)

echo [2/3] build (celeris + celeris_gui, Release)...
"%CMAKE%" --build build-msvc --config Release --target celeris celeris_gui || (echo   BUILD FAILED & exit /b 1)

echo [3/3] physics selftest...
REM MKL/CUDA builds need their DLLs on PATH; CMake writes the dirs here.
set "RTPATH="
if exist "build-msvc\runtime_path.txt" set /p RTPATH=<"build-msvc\runtime_path.txt"
if defined RTPATH set "PATH=%RTPATH%;%PATH%"
echo ----------------------------------------------------------------------
build-msvc\Release\celeris.exe selftest
echo ----------------------------------------------------------------------
echo.
echo CHECK COMPLETE: build OK. Review the selftest numbers above
echo   (expect energy=1.000000, 2D vs grcwa cross-check ^|d^|^<1e-2, Strehl ~0.55).
endlocal
