cmake_minimum_required(VERSION 3.28)
project(CELERIS CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
)
target_include_directories(celeris_core PUBLIC src)
target_link_libraries(celeris_core PUBLIC Eigen3::Eigen)
# ---------------------------------------------------------------------------

# --- celeris: demo / validation app -----------------------------------------
add_executable(celeris main.cpp)
target_link_libraries(celeris PRIVATE celeris_core)

# MinGW GCC: std::print/println terminal output lives in libstdc++exp
if(MINGW)
    target_link_libraries(celeris PRIVATE stdc++exp)
endif()
