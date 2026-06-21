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
set(CELERIS_SUMMARY "")  # appended to as features resolve; printed at the end
set(CELERIS_DLLS "")     # MKL/CUDA runtime DLLs to bundle next to the exe

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
    src/celeris/design/phase_profile.cpp
    src/celeris/design/metalens.cpp
    src/celeris/design/optimize.cpp
    src/celeris/design/system_opt.cpp
    src/celeris/design/polar_metalens.cpp
    src/celeris/design/pb_metalens.cpp
    src/celeris/design/achromatic.cpp
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
    target_include_directories(celeris_core PRIVATE "${CELERIS_CUDA_DIR}/include")

    # Far-field propagation KERNEL -- the real GPU win (600-700x). We COMPILE it by
    # invoking nvcc DIRECTLY (a custom command), NOT via enable_language(CUDA):
    # that avoids enable_language's hard configure error and the MSBuild "CUDA
    # <ver>.targets" integration at build time. The kernel is one self-contained TU
    # (no -rdc) and links the STATIC CUDA runtime (cudart_static), so the GPU build
    # adds ZERO CUDA DLLs -- nothing to bundle, nothing on PATH.
    #
    # GATE: we DELIBERATELY do NOT use check_language(CUDA). Under the VS generator
    # that GENERATES AND BUILDS a throwaway CUDA probe project (CheckCUDA.sln /
    # CompilerIdCUDA.vcxproj) -- and on this setup building that probe pops a Windows
    # "How do you want to open this .vcxproj?" shell-association dialog on every
    # (re)configure (ZERO_CHECK reruns configure at build time, so it recurs per
    # build -- annoying enough to burn into an OLED), and the probe is also flaky
    # (intermittent NOTFOUND -> silent CPU fallback). Instead, probe nvcc the cheap,
    # deterministic, GUI-free way: run `nvcc --version`. It confirms the compiler is
    # invokable WITHOUT building any project. nvcc --version needs no host-compiler
    # env, so it never false-negatives on a bare (no-vcvars) shell; the kernel TU is
    # compiled at BUILD time by the custom command below, inside MSBuild's VC env.
    find_program(CELERIS_NVCC NAMES nvcc HINTS "${CELERIS_CUDA_DIR}/bin")
    set(_cel_nvcc_ok OFF)
    if(CELERIS_NVCC)
        execute_process(COMMAND "${CELERIS_NVCC}" --version
                        RESULT_VARIABLE _cel_nvcc_rc OUTPUT_QUIET ERROR_QUIET)
        if(_cel_nvcc_rc EQUAL 0)
            set(_cel_nvcc_ok ON)
        endif()
    endif()
    if(_cel_nvcc_ok)
        get_filename_component(_cel_cl_dir "${CMAKE_CXX_COMPILER}" DIRECTORY)
        set(_cel_kobj "${CMAKE_BINARY_DIR}/celeris_propagate_cuda${CMAKE_CXX_OUTPUT_EXTENSION}")
        add_custom_command(
            OUTPUT "${_cel_kobj}"
            COMMAND "${CELERIS_NVCC}" -c
                    "${CMAKE_CURRENT_SOURCE_DIR}/src/celeris/cuda/propagate.cu"
                    -o "${_cel_kobj}" -std=c++17 -O3 --cudart static --use-local-env
                    -allow-unsupported-compiler  # tolerate a cl newer than the toolkit
                    -ccbin "${_cel_cl_dir}"
                    -gencode=arch=compute_75,code=sm_75
                    -gencode=arch=compute_86,code=sm_86
                    -gencode=arch=compute_89,code=sm_89
                    -I "${CMAKE_CURRENT_SOURCE_DIR}/src" -I "${CELERIS_CUDA_DIR}/include"
                    -Xcompiler=/MD,/EHsc,/O2,/DNDEBUG
            DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/celeris/cuda/propagate.cu"
            COMMENT "Compiling CUDA far-field kernel via nvcc (generator-agnostic)"
            VERBATIM)
        set_source_files_properties("${_cel_kobj}" PROPERTIES
            EXTERNAL_OBJECT TRUE GENERATED TRUE)
        target_sources(celeris_core PRIVATE "${_cel_kobj}")
        target_link_libraries(celeris_core PUBLIC "${CELERIS_CUDA_DIR}/lib/x64/cudart_static.lib")
        target_compile_definitions(celeris_core PUBLIC CELERIS_USE_CUDA_KERNELS)
        string(APPEND CELERIS_SUMMARY
            "  GPU backend   : CUDA far-field kernel (direct nvcc, static cudart, any generator)\n")
    else()
        string(APPEND CELERIS_SUMMARY
            "  GPU backend   : CUDA toolkit found but this generator/env can't build it\n"
            "                  (CPU propagation; use scripts/build-cuda.bat for the GPU kernel)\n")
    endif()

    # cuSOLVER GPU eigensolve: OPT-IN (default OFF). It is a recorded DEAD END (no
    # faster than the CPU path) and it drags in cusolver + cublasLt as a load-time
    # dependency (~630 MB of DLLs). Only the `gpubench` command uses it.
    option(CELERIS_USE_CUDA_EIGENSOLVE "Build the (dead-end) cuSOLVER GPU eigensolve + gpubench" OFF)
    if(CELERIS_USE_CUDA_EIGENSOLVE)
        target_sources(celeris_core PRIVATE src/celeris/cuda/eigensolve.cpp)
        target_compile_definitions(celeris_core PUBLIC CELERIS_USE_CUDA)
        target_link_libraries(celeris_core PUBLIC
            "${CELERIS_CUDA_DIR}/lib/x64/cusolver.lib"
            "${CELERIS_CUDA_DIR}/lib/x64/cudart.lib")
        file(GLOB _cuda_dll
            "${CELERIS_CUDA_DIR}/bin/x64/cudart64*.dll"
            "${CELERIS_CUDA_DIR}/bin/x64/cusolver64*.dll"
            "${CELERIS_CUDA_DIR}/bin/x64/cublas64*.dll"
            "${CELERIS_CUDA_DIR}/bin/x64/cublasLt64*.dll")
        list(APPEND CELERIS_DLLS ${_cuda_dll})
        string(APPEND CELERIS_SUMMARY "  GPU eigensolve: cuSOLVER (opt-in; bundles ~630 MB of DLLs)\n")
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
    # Bundle the MKL runtime next to the exe (no PATH at runtime). mkl_rt loads,
    # at runtime, core + the threading layer (intel_thread + libiomp5) + the ISA
    # kernels it dispatches to: both the BLAS/LAPACK kernels (mkl_avx2/avx512/def)
    # AND the vector-math kernels (mkl_vml_*). Miss the vml ones and mkl_rt aborts
    # with "cannot load mkl_vml_avx2.dll". Glob the families so it's version-proof.
    file(GLOB _mkl_dll
        "${CELERIS_MKL_BINDIR}/mkl_rt.*.dll"
        "${CELERIS_MKL_BINDIR}/mkl_core.*.dll"
        "${CELERIS_MKL_BINDIR}/mkl_intel_thread.*.dll"
        "${CELERIS_MKL_BINDIR}/mkl_def.*.dll"
        "${CELERIS_MKL_BINDIR}/mkl_avx*.dll"       # avx2 / avx512 BLAS kernels
        "${CELERIS_MKL_BINDIR}/mkl_vml_*.dll"      # vector-math kernels (all ISAs)
        "${CELERIS_MKL_BINDIR}/mkl_mc*.dll"        # pre-AVX fallback kernels
        "${CELERIS_MKL_BINDIR}/libiomp5md.dll")
    list(APPEND CELERIS_DLLS ${_mkl_dll})
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

# --- Python bindings (pybind11) ---------------------------------------------
# Opt-in: needs a Python with development headers. Enable with
# -DCELERIS_BUILD_PYTHON=ON, or let `pip install .` (scikit-build-core) drive it.
# Builds the `_celeris` extension over celeris_core and STAGES it next to the
# pure-Python `celeris` package under build/python/celeris -- so
# `PYTHONPATH=build/python python -c "import celeris"` works straight from the
# build tree, no install step. The same MKL/CUDA runtime DLLs the exe bundles are
# copied next to the .pyd (Windows resolves a module's deps from its own dir).
option(CELERIS_BUILD_PYTHON "Build the Python bindings (pybind11)" OFF)
if(CELERIS_BUILD_PYTHON)
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
    FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v3.0.4 GIT_SHALLOW ON)
    FetchContent_MakeAvailable(pybind11)

    set(_cel_pydir "${CMAKE_BINARY_DIR}/python/celeris")
    pybind11_add_module(_celeris bindings/python/celeris_module.cpp)
    target_link_libraries(_celeris PRIVATE celeris_core)
    # Pin the extension's output dir across every config (VS multi-config would
    # otherwise scatter it into Release/ Debug/ subfolders).
    set_target_properties(_celeris PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${_cel_pydir}")
    foreach(_cfg DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
        set_target_properties(_celeris PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY_${_cfg} "${_cel_pydir}")
    endforeach()

    # Stage the pure-Python package next to the extension (build-tree usable).
    add_custom_command(TARGET _celeris POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                "${CMAKE_CURRENT_SOURCE_DIR}/python/celeris" "${_cel_pydir}"
        COMMENT "Staging the celeris Python package -> build/python/celeris")
    if(CELERIS_DLLS AND WIN32)
        list(REMOVE_DUPLICATES CELERIS_DLLS)
        add_custom_command(TARGET _celeris POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CELERIS_DLLS} "${_cel_pydir}"
            COMMAND_EXPAND_LISTS VERBATIM
            COMMENT "Bundling MKL/CUDA runtime DLLs next to _celeris")
    endif()

    # scikit-build-core (`pip install .`) sets SKBUILD: install into the wheel.
    if(DEFINED SKBUILD)
        install(TARGETS _celeris LIBRARY DESTINATION celeris)
        install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python/celeris/"
                DESTINATION celeris)
    endif()

    string(APPEND CELERIS_SUMMARY
        "  Python bindings: ON (_celeris -> build/python/celeris)\n")
endif()
# ---------------------------------------------------------------------------

# --- Bundle runtime DLLs next to the exe ------------------------------------
# MKL/CUDA link to DLLs that must be found at load time. Windows searches the
# executable's OWN directory before PATH, so copying them next to celeris.exe /
# celeris_gui.exe makes the build self-contained -- no PATH, no launcher script,
# no per-lab setup. This is the deployment model: ship the exe + its DLLs as one
# folder. (Both exes share a directory, so a single copy onto celeris's dir
# serves both.) copy_if_different => only re-copies when a DLL actually changes.
if(CELERIS_DLLS AND WIN32)
    list(REMOVE_DUPLICATES CELERIS_DLLS)
    add_custom_command(TARGET celeris POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CELERIS_DLLS} "$<TARGET_FILE_DIR:celeris>"
        COMMAND_EXPAND_LISTS VERBATIM
        COMMENT "Bundling MKL/CUDA runtime DLLs next to the executables")
    string(APPEND CELERIS_SUMMARY
        "  Runtime DLLs  : bundled next to the exe (runs with no PATH setup)\n")
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
