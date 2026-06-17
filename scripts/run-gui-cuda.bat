@echo off
REM ============================================================================
REM Launch the GPU-accelerated CELERIS GUI built by scripts\build-cuda.bat.
REM CUDA 13 puts its runtime DLLs (cudart, cusolver, ...) in bin\x64, which must
REM be on PATH for the executable to start. Run from the repo root.
REM ============================================================================
set "CUDADIR=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.3"
set "PATH=%CUDADIR%\bin\x64;%CUDADIR%\bin;%PATH%"
start "" "%~dp0..\build-cuda\celeris_gui.exe"
