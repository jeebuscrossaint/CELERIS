cmake_minimum_required(VERSION 3.28)
project(CELERIS CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Optimization matters enormously here: Eigen is 10-100x slower without -O2 and
# with assertions enabled, and dense eigensolves dominate RCWA runtime. Default
# to Release (-O3 -DNDEBUG) unless the user asks otherwise.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# --- Auto-configuration -----------------------------------------------------
# CELERIS_AUTO (default ON): detect this machine's capabilities and build the
# best variant it can, degrading gracefully -- so a single `cmake -B build` does
# the right thing instead of needing a separate dir per backend. The GPU and MKL
# switches are tri-state: AUTO (detect), ON (require, error if missing), or OFF.
# A build-configuration summary is printed at the end of configure.
option(CELERIS_AUTO "Auto-detect and enable the best available backends" ON)
set(CELERIS_USE_CUDA AUTO CACHE STRING "GPU backend: AUTO, ON, or OFF")
set(CELERIS_USE_MKL  AUTO CACHE STRING "Intel MKL backend: AUTO, ON, or OFF")
set_property(CACHE CELERIS_USE_CUDA PROPERTY STRINGS AUTO ON OFF)
set_property(CACHE CELERIS_USE_MKL  PROPERTY STRINGS AUTO ON OFF)
set(CELERIS_SUMMARY "")       # appended to as features resolve; printed at the end
set(CELERIS_RUNTIME_DIRS "")  # dirs that must be on PATH at runtime (MKL/CUDA DLLs)

# Vectorization: MSVC defaults to the SSE2 x64 baseline, which leaves Eigen's
# dense complex matmuls (the real RCWA hot path -- operator assembly + S-matrix
# recursion, far more costly than the eigensolve) badly under-vectorized. Enable
# AVX2+FMA, a big free speedup on any modern CPU. Opt out with -DCELERIS_AVX2=OFF
# for portability to pre-2013 hardware.
option(CELERIS_AVX2 "Build with AVX2/FMA (fast; needs a Haswell-or-newer CPU)" ON)
if(CELERIS_AVX2)
    # Scope to C++ only: a bare /arch:AVX2 leaks into the nvcc command line and
    # breaks the CUDA kernel compile (nvcc doesn't accept it directly).
    if(MSVC)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/arch:AVX2>)
    else()
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-mavx2>
                            $<$<COMPILE_LANGUAGE:CXX>:-mfma>)
    endif()
    string(APPEND CELERIS_SUMMARY "  Vectorization : AVX2/FMA\n")
else()
    string(APPEND CELERIS_SUMMARY "  Vectorization : SSE2 baseline (AVX2 off)\n")
endif()

# --- External libraries -----------------------------------------------------
include(FetchContent)

# Eigen: header-only linear algebra. Downloaded at configure time, no install.
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    ON
)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_BUILD_TESTING OFF)
set(BUILD_TESTING OFF)
FetchContent_MakeAvailable(Eigen)
# ---------------------------------------------------------------------------

# --- celeris_core: the GUI-agnostic engine library --------------------------
# Everything physics/numerics lives here and returns plain data. Any front end
# (CLI, native GUI, Python bindings) links this and never the other way round.
add_library(celeris_core
    src/celeris/materials/material.cpp
    src/celeris/materials/database.cpp
    src/celeris/optics/tmm.cpp
    src/celeris/rcwa/grating1d.cpp
    src/celeris/rcwa/rcwa1d.cpp
    src/celeris/rcwa/rcwa_stack1d.cpp
    src/celeris/rcwa/rcwa2d.cpp
    src/celeris/rcwa/eig.cpp
    src/celeris/design/metalens.cpp
    src/celeris/design/optimize.cpp
    src/celeris/design/system_opt.cpp
    src/celeris/design/polar_metalens.cpp
    src/celeris/analysis/focal.cpp
    src/celeris/analysis/chromatic.cpp
    src/celeris/analysis/tolerance.cpp
    src/celeris/analysis/field.cpp
    src/celeris/analysis/wavefront.cpp
    src/celeris/analysis/mtf.cpp
    src/celeris/analysis/throughfocus.cpp
    src/celeris/analysis/polarization.cpp
    src/celeris/io/gds.cpp
    src/celeris/io/image.cpp
    src/celeris/io/material_csv.cpp
)
find_package(Threads REQUIRED)
target_include_directories(celeris_core PUBLIC src)
target_link_libraries(celeris_core PUBLIC Eigen3::Eigen Threads::Threads)

# --- CUDA GPU backend (cuSOLVER eigensolve + far-field kernel) ---------------
# AUTO: detect a CUDA toolkit (CUDA_PATH / CUDA_HOME env, or scoop) and enable if
# present. ON: require it. OFF: skip. The host cuSOLVER path builds with any
# generator; the far-field *kernel* (the 600-700x win) needs nvcc, which the VS
# generator can't drive without the VS CUDA integration -- so on a plain VS build
# you get a clear "use the Ninja path (scripts/build-cuda.bat)" note.
set(_cel_cuda_dir "")
foreach(_cand "$ENV{CUDA_PATH}" "$ENV{CUDA_HOME}" "$ENV{USERPROFILE}/scoop/apps/cuda/current")
    if(_cand AND EXISTS "${_cand}/lib/x64/cudart.lib" AND NOT _cel_cuda_dir)
        set(_cel_cuda_dir "${_cand}")
    endif()
endforeach()

set(_cel_cuda OFF)
if(CELERIS_USE_CUDA STREQUAL "OFF")
    string(APPEND CELERIS_SUMMARY "  GPU backend   : off (CELERIS_USE_CUDA=OFF)\n")
elseif(_cel_cuda_dir)
    set(_cel_cuda ON)
elseif(CELERIS_USE_CUDA STREQUAL "ON")
    message(FATAL_ERROR "CELERIS_USE_CUDA=ON but no CUDA toolkit found "
                        "(set CUDA_PATH or install the toolkit).")
else()
    string(APPEND CELERIS_SUMMARY "  GPU backend   : CPU only (no CUDA toolkit found)\n")
endif()

if(_cel_cuda)
    set(CELERIS_CUDA_DIR "${_cel_cuda_dir}"
        CACHE PATH "CUDA toolkit directory (contains include/ and lib/x64/)" FORCE)
    list(APPEND CELERIS_RUNTIME_DIRS "${CELERIS_CUDA_DIR}/bin")  # cudart/cusolver DLLs
    # The GPU eigensolve is pure HOST code (cuSOLVER host API), so we compile it
    # with the normal C++ compiler and just link the CUDA import libs. No nvcc,
    # no CUDA language, no CCCL headers (which scoop's CUDA package omits), and
    # it builds with the plain Visual Studio generator.
    target_sources(celeris_core PRIVATE src/celeris/cuda/eigensolve.cpp)
    target_compile_definitions(celeris_core PUBLIC CELERIS_USE_CUDA)
    target_include_directories(celeris_core PRIVATE "${CELERIS_CUDA_DIR}/include")
    target_link_libraries(celeris_core PUBLIC
        "${CELERIS_CUDA_DIR}/lib/x64/cusolver.lib"
        "${CELERIS_CUDA_DIR}/lib/x64/cudart.lib")

    # Device kernels (the far-field propagation) need nvcc. Build them in a
    # SEPARATE CUDA-only static library so they don't inherit the project's
    # C++23 dialect (nvcc would reject "CUDA23"); the kernel only needs C++17.
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES 75 86 89)  # Turing/Ampere/Ada
        endif()
        add_library(celeris_cuda STATIC src/celeris/cuda/propagate.cu)
        target_include_directories(celeris_cuda PRIVATE src "${CELERIS_CUDA_DIR}/include")
        set_target_properties(celeris_cuda PROPERTIES
            CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON
            CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
        target_link_libraries(celeris_core PUBLIC celeris_cuda)
        target_compile_definitions(celeris_core PUBLIC CELERIS_USE_CUDA_KERNELS)
        string(APPEND CELERIS_SUMMARY
            "  GPU backend   : CUDA cuSOLVER + far-field kernel (nvcc OK)\n")
    else()
        string(APPEND CELERIS_SUMMARY
            "  GPU backend   : CUDA cuSOLVER only -- kernel needs nvcc; this\n"
            "                  generator can't drive it, use scripts/build-cuda.bat (Ninja)\n")
    endif()
endif()
# ---------------------------------------------------------------------------

# --- Optional LAPACK eigensolve backend -------------------------------------
# The RCWA per-layer eigenproblem (general complex zgeev) dominates the library
# build. eig.cpp can dispatch it to LAPACK instead of Eigen's header-only solver,
# statically linked so the exe stays a single self-contained file.
#
# MEASURED (this machine): the only Fortran-free LAPACK available here (CLAPACK,
# f2c reference LAPACK, even on top of OpenBLAS's fast BLAS) is ~1.5x SLOWER than
# Eigen -- the bottleneck is the *reference* zgeev algorithm, not the BLAS. A
# real speedup needs a re-optimized LAPACK, i.e. Intel MKL. So this flag is a
# plug-in point for MKL (or any optimized LAPACK exposing zgeev_), NOT a win with
# reference LAPACK. Default OFF -> Eigen, which is the fastest option here.
option(CELERIS_USE_OPENBLAS "Link an external LAPACK (e.g. MKL) for zgeev; needs an OPTIMIZED LAPACK to beat Eigen" OFF)
if(CELERIS_USE_OPENBLAS)
    find_library(CEL_LAPACK NAMES lapack mkl_rt REQUIRED)
    find_library(CEL_BLAS NAMES openblas blas mkl_rt)
    find_library(CEL_F2C NAMES libf2c f2c)
    target_link_libraries(celeris_core PUBLIC ${CEL_LAPACK})
    if(CEL_BLAS)
        target_link_libraries(celeris_core PUBLIC ${CEL_BLAS})
    endif()
    if(CEL_F2C)
        target_link_libraries(celeris_core PUBLIC ${CEL_F2C})
    endif()
    target_compile_definitions(celeris_core PUBLIC CELERIS_USE_LAPACK)
    message(STATUS "CELERIS: external LAPACK zgeev enabled -- ${CEL_LAPACK}")
endif()

# Intel MKL turbo. MEASURED: the RCWA cost is dominated NOT by the eigensolve but
# by Eigen's dense complex matmuls / inverses (operator assembly + S-matrix
# recursion). Routing ALL of Eigen's dense BLAS/LAPACK through MKL
# (EIGEN_USE_MKL_ALL) plus the zgeev shim is ~11x over the SSE2 baseline at high M.
# AUTO: search the common spots (a `pip install mkl-devel` python prefix's
# Library/, $MKLROOT, $CONDA_PREFIX/Library, oneAPI) and enable if found; with
# CELERIS_FETCH_MKL=ON, pip-install it first if missing. Enabling MKL means the
# mkl_rt DLL must be on PATH at runtime, so the dependency-free AVX2 build stays
# the fallback. Skipped if the (legacy) explicit OPENBLAS LAPACK is in use.
option(CELERIS_FETCH_MKL "If MKL is not found, pip-install mkl-devel/mkl-include" OFF)
set(_cel_mkl OFF)
if(NOT CELERIS_USE_OPENBLAS AND NOT CELERIS_USE_MKL STREQUAL "OFF")
    # Hint search at conda, oneAPI, MKLROOT, and any Python (pip mkl-devel) prefix.
    # FindPython3 locates the interpreter far more reliably than find_program
    # (registry, standard install dirs) -- needed both to read its prefix and to
    # drive CELERIS_FETCH_MKL.
    find_package(Python3 QUIET COMPONENTS Interpreter)
    set(_cel_py "${Python3_EXECUTABLE}")
    set(_mkl_hints
        "$ENV{MKLROOT}"
        "$ENV{CONDA_PREFIX}/Library"
        "$ENV{ONEAPI_ROOT}/mkl/latest"
        "C:/Program Files (x86)/Intel/oneAPI/mkl/latest"
        "$ENV{USERPROFILE}/scoop/apps/python/current/Library")  # scoop python pip-mkl
    if(_cel_py)
        execute_process(COMMAND "${_cel_py}" -c "import sys;print(sys.prefix)"
            OUTPUT_VARIABLE _cel_pyprefix OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        if(_cel_pyprefix)
            list(APPEND _mkl_hints "${_cel_pyprefix}/Library" "${_cel_pyprefix}")
        endif()
    endif()
    find_path(MKL_INC mkl.h HINTS ${_mkl_hints} PATH_SUFFIXES include)
    find_library(MKL_RT NAMES mkl_rt HINTS ${_mkl_hints} PATH_SUFFIXES lib)
    if(NOT (MKL_INC AND MKL_RT) AND CELERIS_FETCH_MKL AND _cel_py)
        message(STATUS "CELERIS: MKL not found -- fetching via pip (CELERIS_FETCH_MKL=ON)...")
        execute_process(COMMAND "${_cel_py}" -m pip install mkl-devel mkl-include)
        find_path(MKL_INC mkl.h HINTS ${_mkl_hints} PATH_SUFFIXES include)
        find_library(MKL_RT NAMES mkl_rt HINTS ${_mkl_hints} PATH_SUFFIXES lib)
    endif()
    if(MKL_INC AND MKL_RT)
        set(_cel_mkl ON)
    elseif(CELERIS_USE_MKL STREQUAL "ON")
        message(FATAL_ERROR "CELERIS_USE_MKL=ON but MKL not found. Try "
                            "`pip install mkl-devel mkl-include` or -DCELERIS_FETCH_MKL=ON.")
    endif()
endif()

if(_cel_mkl)
    target_include_directories(celeris_core PUBLIC ${MKL_INC})
    target_link_libraries(celeris_core PUBLIC ${MKL_RT})
    target_compile_definitions(celeris_core PUBLIC CELERIS_USE_LAPACK EIGEN_USE_MKL_ALL)
    # Runtime DLL lives in the sibling bin/ of the import lib's lib/ dir.
    get_filename_component(_mkl_libdir "${MKL_RT}" DIRECTORY)
    get_filename_component(_mkl_root "${_mkl_libdir}" DIRECTORY)
    set(CELERIS_MKL_BINDIR "${_mkl_root}/bin")
    list(APPEND CELERIS_RUNTIME_DIRS "${CELERIS_MKL_BINDIR}")  # mkl_rt + kernels
    string(APPEND CELERIS_SUMMARY "  Eigensolve+BLAS: Intel MKL (multithreaded)\n")
elseif(CELERIS_USE_OPENBLAS)
    string(APPEND CELERIS_SUMMARY "  Eigensolve+BLAS: external LAPACK (CELERIS_USE_OPENBLAS)\n")
else()
    string(APPEND CELERIS_SUMMARY "  Eigensolve+BLAS: Eigen (header-only)\n")
    string(APPEND CELERIS_SUMMARY "                  for ~3x more: pip install mkl-devel mkl-include (or -DCELERIS_FETCH_MKL=ON)\n")
endif()
# ---------------------------------------------------------------------------

# --- celeris: demo / validation app -----------------------------------------
add_executable(celeris main.cpp)
target_link_libraries(celeris PRIVATE celeris_core)

# MinGW GCC: std::print/println terminal output lives in libstdc++exp
if(MINGW)
    target_link_libraries(celeris PRIVATE stdc++exp)
endif()

# --- Desktop GUI (Dear ImGui + GLFW + OpenGL) --------------------------------
# On by default -- CELERIS is a desktop design tool, so the GUI is part of "the
# best build". Dependencies (GLFW + Dear ImGui) are fetched automatically (MIT,
# no install). Disable for CLI-only / headless / CI builds with
# -DCELERIS_BUILD_GUI=OFF (headless Linux would also need X11/OpenGL dev libs).
option(CELERIS_BUILD_GUI "Build the Dear ImGui desktop GUI" ON)
if(CELERIS_BUILD_GUI)
    include(FetchContent)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(glfw GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG 3.4 GIT_SHALLOW ON)
    FetchContent_MakeAvailable(glfw)

    # Dear ImGui ships no CMakeLists — fetch sources and build them ourselves.
    # Use the docking branch for a dockable/tabbable multi-window workspace.
    FetchContent_Declare(imgui GIT_REPOSITORY https://github.com/ocornut/imgui.git
        GIT_TAG docking GIT_SHALLOW ON)
    FetchContent_MakeAvailable(imgui)

    find_package(OpenGL REQUIRED)
    add_library(celeris_imgui STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp)
    target_include_directories(celeris_imgui PUBLIC
        ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)
    target_link_libraries(celeris_imgui PUBLIC glfw OpenGL::GL)

    add_executable(celeris_gui
        gui/main_gui.cpp
        gui/theme.cpp
        gui/textures.cpp
        gui/app_state.cpp
        gui/workers.cpp)
    target_link_libraries(celeris_gui PRIVATE celeris_core celeris_imgui)

    # Short git commit -> window title (captured at configure time).
    execute_process(COMMAND git rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE CEL_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(NOT CEL_GIT_COMMIT)
        set(CEL_GIT_COMMIT "dev")
    endif()
    target_compile_definitions(celeris_gui PRIVATE
        CELERIS_GIT_COMMIT="${CEL_GIT_COMMIT}")

    if(WIN32)
        # psapi: process CPU/mem monitor; dwmapi: dark-mode native title bar.
        target_link_libraries(celeris_gui PRIVATE psapi dwmapi)
    endif()
    if(MINGW)
        target_link_libraries(celeris_gui PRIVATE stdc++exp)
    endif()
    string(APPEND CELERIS_SUMMARY "  Desktop GUI   : ON (celeris_gui)\n")
else()
    string(APPEND CELERIS_SUMMARY "  Desktop GUI   : off (enable with -DCELERIS_BUILD_GUI=ON)\n")
endif()

# --- Runtime DLL path (MKL / CUDA) ------------------------------------------
# MKL and CUDA link to DLLs that must be on PATH at runtime. Rather than copy
# ~300 MB of MKL DLLs into every build dir, emit (a) runtime_path.txt -- the
# semicolon-joined dirs, read by scripts/check.bat -- and (b) a celeris-run.bat
# launcher that prepends them to PATH and forwards to celeris.exe, so the exe
# "just works" without the user managing PATH.
if(CELERIS_RUNTIME_DIRS)
    list(REMOVE_DUPLICATES CELERIS_RUNTIME_DIRS)
    set(_rt_join "${CELERIS_RUNTIME_DIRS}")  # quoted list -> "dir1;dir2"
    file(WRITE "${CMAKE_BINARY_DIR}/runtime_path.txt" "${_rt_join}")
    if(WIN32)
        string(REPLACE "/" "\\" _rt_win "${_rt_join}")
        file(WRITE "${CMAKE_BINARY_DIR}/celeris-run.bat"
            "@echo off\r\n"
            "set \"PATH=${_rt_win};%PATH%\"\r\n"
            "set \"EXE=%~dp0Release\\celeris.exe\"\r\n"
            "if not exist \"%EXE%\" set \"EXE=%~dp0celeris.exe\"\r\n"
            "\"%EXE%\" %*\r\n")
    endif()
    string(APPEND CELERIS_SUMMARY
        "  Runtime PATH  : auto via celeris-run.bat (or `call`-set from runtime_path.txt)\n")
endif()
# ---------------------------------------------------------------------------

# --- Build configuration summary --------------------------------------------
message(STATUS "")
message(STATUS "================ CELERIS build configuration ================")
# message() collapses the leading spaces of each line, so print the block raw.
string(REGEX REPLACE "\n$" "" _cel_sum "${CELERIS_SUMMARY}")
message("${_cel_sum}")
message(STATUS "=============================================================")
message(STATUS "")
