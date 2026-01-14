cmake_minimum_required(VERSION 3.15)

project(OrbFit
    VERSION 1.0
    DESCRIPTION "High-precision Asteroid Dynamics Library"
    HOMEPAGE_URL "https://github.com/ITALOccult/AstDyn"
    LANGUAGES CXX C
)

# ============================================================================
# Project Options
# ============================================================================

option(ASTDYN_BUILD_TESTS "Build unit tests" ON)
option(ASTDYN_BUILD_EXAMPLES "Build examples" ON)
option(ASTDYN_BUILD_DOCS "Build documentation" OFF)
option(ASTDYN_USE_SPICE "Use SPICE toolkit for ephemerides" OFF)
option(ASTDYN_ENABLE_PROFILING "Enable profiling support" OFF)
option(ASTDYN_BUILD_SHARED "Build shared library in addition to static" ON)

# ============================================================================
# C++ Standard and Compiler Settings
# ============================================================================

set(CMAKE_CXX_STANDARD 17)
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

# Optional: SPICE Toolkit
if(ASTDYN_USE_SPICE)
    find_package(CSPICE QUIET)
    if(CSPICE_FOUND)
        message(STATUS "CSPICE found: ${CSPICE_INCLUDE_DIRS}")
        add_compile_definitions(ASTDYN_USE_SPICE)
    else()
        message(WARNING "CSPICE not found. Ephemeris functionality will be limited.")
        set(ASTDYN_USE_SPICE OFF)
    endif()
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
message(STATUS "║ Version:          ${PROJECT_VERSION}")
message(STATUS "║ Build Type:       ${CMAKE_BUILD_TYPE}")
message(STATUS "║ C++ Standard:     C++${CMAKE_CXX_STANDARD}")
message(STATUS "║ Compiler:         ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "╠════════════════════════════════════════════════════════════╣")
message(STATUS "║ Options:")
message(STATUS "║   Build Tests:    ${ASTDYN_BUILD_TESTS}")
message(STATUS "║   Build Examples: ${ASTDYN_BUILD_EXAMPLES}")
message(STATUS "║   Build Docs:     ${ASTDYN_BUILD_DOCS}")
message(STATUS "║   Use SPICE:      ${ASTDYN_USE_SPICE}")
message(STATUS "║   Profiling:      ${ASTDYN_ENABLE_PROFILING}")
message(STATUS "╠════════════════════════════════════════════════════════════╣")
message(STATUS "║ Dependencies:")
message(STATUS "║   Eigen3:         ${EIGEN3_VERSION_STRING}")
message(STATUS "║   Boost:          ${Boost_VERSION_STRING}")
if(ASTDYN_USE_SPICE AND CSPICE_FOUND)
message(STATUS "║   CSPICE:         Found")
endif()
if(ASTDYN_BUILD_TESTS)
message(STATUS "║   Google Test:    Enabled")
endif()
message(STATUS "╚════════════════════════════════════════════════════════════╝")
message(STATUS "")
