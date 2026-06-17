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
    src/celeris/design/metalens.cpp
    src/celeris/design/optimize.cpp
    src/celeris/analysis/focal.cpp
    src/celeris/analysis/chromatic.cpp
    src/celeris/analysis/tolerance.cpp
    src/celeris/io/gds.cpp
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
endif()
# ---------------------------------------------------------------------------

# --- celeris: demo / validation app -----------------------------------------
add_executable(celeris main.cpp)
target_link_libraries(celeris PRIVATE celeris_core)

# MinGW GCC: std::print/println terminal output lives in libstdc++exp
if(MINGW)
    target_link_libraries(celeris PRIVATE stdc++exp)
endif()
