cmake_minimum_required(VERSION 3.10)

# Project name
project(ASTDYN LANGUAGES CXX)

# Require C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen (header‑only) – if not found, fall back to include path
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(nlohmann_json 3.2.0 REQUIRED)
# Find GoogleTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Include directories for the AstDyn library
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/astdyn/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/astdyn)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/astdyn/external/ioc_gaialib/include)

# Add source files (you can extend this list as needed)
file(GLOB_RECURSE ASTDYN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/astdyn/src/*.cpp")
list(FILTER ASTDYN_SOURCES EXCLUDE REGEX ".*astdyn_main\\.cpp$")
list(FILTER ASTDYN_SOURCES EXCLUDE REGEX ".*astdyn_convert\\.cpp$")
list(FILTER ASTDYN_SOURCES EXCLUDE REGEX ".*astdyn/src/python/.*")
add_library(astdyn STATIC ${ASTDYN_SOURCES})

target_link_libraries(astdyn Eigen3::Eigen nlohmann_json::nlohmann_json)

# Example executable (propagation demo)
add_executable(propagate_example astdyn/examples/propagate_79148.cpp)
target_link_libraries(propagate_example astdyn Eigen3::Eigen)

# Benchmarks
add_subdirectory(astdyn/benchmarks)

# Enable testing
enable_testing()
