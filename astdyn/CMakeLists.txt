cmake_minimum_required(VERSION 3.15)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

project(OrbFit
    VERSION 0.9.0
    DESCRIPTION "High-precision Asteroid Dynamics Library (beta 0.9)"
    HOMEPAGE_URL "https://github.com/ITALOccult/AstDyn"
    LANGUAGES CXX C
)

# ============================================================================
# Project Options
# ============================================================================

option(ASTDYN_BUILD_TESTS "Build unit tests" OFF)
option(ASTDYN_BUILD_EXAMPLES "Build examples" ON)
option(ASTDYN_BUILD_DOCS "Build documentation" OFF)
option(ASTDYN_USE_SPICE "Use SPICE toolkit for ephemerides" OFF)
option(ASTDYN_ENABLE_PROFILING "Enable profiling support" OFF)
option(ASTDYN_BUILD_SHARED "Build shared library in addition to static" ON)
option(ASTDYN_BUILD_CATALOG "Build Gaia DR3 stellar catalog module (requires local catalog data)" OFF)
option(ASTDYN_USE_OPENMP "Enable OpenMP parallelization" ON)

# ============================================================================
# C++ Standard and Compiler Settings
# ============================================================================

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile commands for IDE integration
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ============================================================================
# Build Type Configuration
# ============================================================================

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ============================================================================
# Compiler Flags
# ============================================================================

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        -Wall 
        -Wextra 
        -Wpedantic
        -Wno-unused-parameter
    )
    
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g -O0)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(-O3 -march=native)
    endif()
    
    if(ASTDYN_ENABLE_PROFILING)
        add_compile_options(-pg)
        add_link_options(-pg)
    endif()
    
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options(/W4 /WX)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(/O2)
    endif()
endif()

# ============================================================================
# Dependencies
# ============================================================================

# Include CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Find or fetch Eigen3
find_package(Eigen3 3.4 QUIET)
if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found, fetching from GitLab...")
    include(FetchContent)
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(Eigen3)
endif()

# Find Boost (filesystem, program_options, date_time)
find_package(Boost 1.70 REQUIRED COMPONENTS 
    filesystem 
    program_options
    date_time
)

# Find CURL for networking
find_package(CURL REQUIRED)

# Find or fetch nlohmann/json
find_package(nlohmann_json QUIET)
if(NOT nlohmann_json_FOUND)
    message(STATUS "nlohmann_json not found, fetching...")
    include(FetchContent)
    FetchContent_Declare(
        json
        URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
    )
    FetchContent_MakeAvailable(json)
endif()

# Find OpenMP
if(ASTDYN_USE_OPENMP)
    if(APPLE)
        # Search for Homebrew libomp
        if(EXISTS "/opt/homebrew/opt/libomp")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor;-fopenmp")
            set(OpenMP_CXX_INCLUDE_DIRS "/opt/homebrew/opt/libomp/include")
            set(OpenMP_CXX_LIB_NAMES "omp")
            set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
            set(OpenMP_FOUND TRUE)
            message(STATUS "OpenMP found via Homebrew: /opt/homebrew/opt/libomp")
        endif()
    endif()
    
    if(NOT OpenMP_FOUND)
        find_package(OpenMP)
    endif()
    
    if(OpenMP_FOUND)
        message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
        if(NOT TARGET OpenMP::OpenMP_CXX)
            add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
            set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
                INTERFACE_COMPILE_OPTIONS "${OpenMP_CXX_FLAGS}"
                INTERFACE_INCLUDE_DIRECTORIES "${OpenMP_CXX_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${OpenMP_omp_LIBRARY}"
            )
        endif()
    else()
        message(WARNING "OpenMP not found, parallelization will be disabled")
    endif()
endif()

# Python Bindings (pybind11)
option(ASTDYN_BUILD_PYTHON "Build Python bindings" OFF)
if(ASTDYN_BUILD_PYTHON)
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG v2.11.1
    )
    FetchContent_MakeAvailable(pybind11)
endif()

# Catalog bridge symbols are referenced by core astrometry code even when the
# local Gaia dataset is not enabled at runtime.
# ioc_gaialib has been internalized into src/catalog/gaia/ and include/astdyn/catalog/gaia/
find_package(SQLite3 REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(ZLIB REQUIRED)

# Optional: Profiling
if(ASTDYN_ENABLE_PROFILING)
    add_compile_definitions(ASTDYN_ENABLE_PROFILING)
endif()

# Google Test (only if building tests)
if(ASTDYN_BUILD_TESTS)
    find_package(GTest QUIET)
    if(NOT GTest_FOUND)
        message(STATUS "Google Test not found, fetching from GitHub...")
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.14.0
            GIT_SHALLOW TRUE
        )
        # For Windows: Prevent overriding parent project's compiler/linker settings
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    endif()
    enable_testing()
    include(GoogleTest)
endif()

# ============================================================================
# Project Configuration
# ============================================================================

# Generate version header
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Version.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/include/astdyn/Version.hpp"
    @ONLY
)

# Generate config header
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/include/astdyn/Config.hpp"
    @ONLY
)

# ============================================================================
# Include Directories
# ============================================================================

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/include
)

# ============================================================================
# Subdirectories
# ============================================================================

add_subdirectory(src)

if(ASTDYN_BUILD_TESTS)
    add_subdirectory(tests)
endif()

if(ASTDYN_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# add_subdirectory(benchmarks)
add_subdirectory(tools)

if(ASTDYN_BUILD_PYTHON)
    add_subdirectory(src/python)
endif()

if(ASTDYN_BUILD_DOCS)
    add_subdirectory(docs)
endif()

# ============================================================================
# Installation
# ============================================================================

include(GNUInstallDirs)

# Install headers
install(DIRECTORY include/astdyn
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)

# Install generated headers
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/astdyn
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)

# Install CMake package config
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/AstDynConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/AstDynConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/OrbFit
)

# ============================================================================
# Summary
# ============================================================================

message(STATUS "")
message(STATUS "╔════════════════════════════════════════════════════════════╗")
message(STATUS "║              OrbFit C++ Configuration Summary              ║")
message(STATUS "╠════════════════════════════════════════════════════════════╣")
message(STATUS "║ Version:          beta 0.9")
message(STATUS "║ Build Type:       ${CMAKE_BUILD_TYPE}")
message(STATUS "║ C++ Standard:     C++${CMAKE_CXX_STANDARD}")
message(STATUS "║ Compiler:         ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "╠════════════════════════════════════════════════════════════╣")
message(STATUS "║ Options:")
message(STATUS "║   Build Tests:    ${ASTDYN_BUILD_TESTS}")
message(STATUS "║   Build Examples: ${ASTDYN_BUILD_EXAMPLES}")
message(STATUS "║   Build Docs:     ${ASTDYN_BUILD_DOCS}")
message(STATUS "║   Profiling:      ${ASTDYN_ENABLE_PROFILING}")
message(STATUS "║   Gaia Catalog:   ${ASTDYN_BUILD_CATALOG}")
message(STATUS "╠════════════════════════════════════════════════════════════╣")
message(STATUS "║ Dependencies:")
message(STATUS "║   Eigen3:         ${EIGEN3_VERSION_STRING}")
message(STATUS "║   Boost:          ${Boost_VERSION_STRING}")
if(ASTDYN_BUILD_TESTS)
message(STATUS "║   Google Test:    Enabled")
endif()
message(STATUS "╚════════════════════════════════════════════════════════════╝")
message(STATUS "")
