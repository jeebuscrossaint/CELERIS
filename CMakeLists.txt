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

# Vectorization: MSVC defaults to the SSE2 x64 baseline, which leaves Eigen's
# dense complex matmuls (the real RCWA hot path -- operator assembly + S-matrix
# recursion, far more costly than the eigensolve) badly under-vectorized. Enable
# AVX2+FMA, a big free speedup on any modern CPU. Opt out with -DCELERIS_AVX2=OFF
# for portability to pre-2013 hardware.
option(CELERIS_AVX2 "Build with AVX2/FMA (fast; needs a Haswell-or-newer CPU)" ON)
if(CELERIS_AVX2)
    if(MSVC)
        add_compile_options(/arch:AVX2)
    else()
        add_compile_options(-mavx2 -mfma)
    endif()
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

# --- Optional CUDA GPU backend (eigensolve via cuSOLVER) ---------------------
# Enable with -DCELERIS_USE_CUDA=ON. On Windows this needs the MSVC toolchain
# plus the VS generator toolset:
#   -G "Visual Studio 17 2022" -A x64
#   -DCMAKE_GENERATOR_TOOLSET="cuda=<scoop>/apps/cuda/current"
option(CELERIS_USE_CUDA "Build the CUDA GPU eigensolve backend" OFF)
if(CELERIS_USE_CUDA)
    # The GPU eigensolve is pure HOST code (cuSOLVER host API), so we compile it
    # with the normal C++ compiler and just link the CUDA import libs. No nvcc,
    # no CUDA language, no CCCL headers (which scoop's CUDA package omits), and
    # it builds with the plain Visual Studio generator.
    set(CELERIS_CUDA_DIR "$ENV{USERPROFILE}/scoop/apps/cuda/current"
        CACHE PATH "CUDA toolkit directory (contains include/ and lib/x64/)")
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
        message(STATUS "CELERIS: CUDA device kernels enabled (nvcc found)")
    else()
        message(STATUS "CELERIS: nvcc not found - GPU propagation kernel disabled "
                       "(host-only cuSOLVER path still built)")
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

# Intel MKL turbo (opt-in). MEASURED: the RCWA cost is dominated NOT by the
# eigensolve but by Eigen's dense complex matmuls / inverses (operator assembly +
# S-matrix recursion). So MKL-zgeev alone buys only ~15%; the real win is routing
# ALL of Eigen's dense BLAS/LAPACK through MKL (EIGEN_USE_MKL_ALL) PLUS the zgeev
# shim. Needs MKL on the include/lib path (e.g. `pip install mkl-devel mkl-include`
# gives <prefix>/Library/{include,lib,bin}); the mkl_rt DLL must be on PATH at run.
# Trades the single-self-contained-binary property for speed -- hence opt-in; the
# default build is AVX2 Eigen (dependency-free, already ~2.7x over the SSE2 base).
option(CELERIS_USE_MKL "Route Eigen BLAS/LAPACK + zgeev through Intel MKL (needs MKL on PATH)" OFF)
if(CELERIS_USE_MKL)
    find_path(MKL_INC mkl.h REQUIRED)
    find_library(MKL_RT NAMES mkl_rt REQUIRED)
    target_include_directories(celeris_core PUBLIC ${MKL_INC})
    target_link_libraries(celeris_core PUBLIC ${MKL_RT})
    target_compile_definitions(celeris_core PUBLIC CELERIS_USE_LAPACK EIGEN_USE_MKL_ALL)
    message(STATUS "CELERIS: Intel MKL backend enabled -- ${MKL_RT} (inc ${MKL_INC})")
endif()
# ---------------------------------------------------------------------------

# --- celeris: demo / validation app -----------------------------------------
add_executable(celeris main.cpp)
target_link_libraries(celeris PRIVATE celeris_core)

# MinGW GCC: std::print/println terminal output lives in libstdc++exp
if(MINGW)
    target_link_libraries(celeris PRIVATE stdc++exp)
endif()

# --- Optional desktop GUI (Dear ImGui + GLFW + OpenGL) -----------------------
# Enable with -DCELERIS_BUILD_GUI=ON. Dependencies fetched automatically (MIT
# licensed, no install). Builds the `celeris_gui` interactive app.
option(CELERIS_BUILD_GUI "Build the Dear ImGui desktop GUI" OFF)
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
endif()
