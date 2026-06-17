@echo off
REM ============================================================================
REM Build CELERIS with the CUDA GPU backend on Windows (MSVC + Ninja + scoop).
REM
REM We set the MSVC + Windows SDK environment EXPLICITLY rather than calling
REM vcvars64.bat: on this machine vcvars fails to add the Windows SDK to LIB/
REM PATH (kernel32.lib / rc.exe / mt.exe missing), which breaks Ninja builds.
REM Setting the paths directly is fully deterministic.
REM
REM If you upgrade VS or the Windows SDK, update VCVER / SDKVER below.
REM Run from the repo root:  scripts\build-cuda.bat
REM ============================================================================
setlocal

set "VSROOT=C:\Program Files\Microsoft Visual Studio\2022\Community"
set "VCVER=14.44.35207"
set "SDKROOT=C:\Program Files (x86)\Windows Kits\10"
set "SDKVER=10.0.26100.0"

REM Full official CUDA Toolkit (has CCCL + nvcc; scoop's package omits CCCL).
set "CUDADIR=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.3"

set "VCT=%VSROOT%\VC\Tools\MSVC\%VCVER%"
set "INCLUDE=%VCT%\include;%SDKROOT%\Include\%SDKVER%\ucrt;%SDKROOT%\Include\%SDKVER%\um;%SDKROOT%\Include\%SDKVER%\shared;%SDKROOT%\Include\%SDKVER%\winrt"
set "LIB=%VCT%\lib\x64;%SDKROOT%\Lib\%SDKVER%\ucrt\x64;%SDKROOT%\Lib\%SDKVER%\um\x64"
set "PATH=%VCT%\bin\Hostx64\x64;%SDKROOT%\bin\%SDKVER%\x64;%CUDADIR%\bin;%PATH%"

set "CMAKE=%USERPROFILE%\scoop\shims\cmake.exe"
set "NINJA=%USERPROFILE%\scoop\shims\ninja.exe"
set "NVCC=%CUDADIR%\bin\nvcc.exe"

"%CMAKE%" -B build-cuda -G Ninja ^
  -DCMAKE_MAKE_PROGRAM="%NINJA%" ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_C_COMPILER=cl ^
  -DCMAKE_CXX_COMPILER=cl ^
  -DCMAKE_CUDA_COMPILER="%NVCC%" ^
  -DCMAKE_CUDA_ARCHITECTURES=89 ^
  -DCELERIS_USE_CUDA=ON ^
  -DCELERIS_CUDA_DIR="%CUDADIR%" || exit /b 1

"%CMAKE%" --build build-cuda --target celeris || exit /b 1
echo BUILD-CUDA-OK
