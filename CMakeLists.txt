cmake_minimum_required(VERSION 3.10)

# Project name
project(ASTDYN LANGUAGES CXX)

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen (header‑only) – if not found, fall back to include path
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(nlohmann_json 3.2.0 REQUIRED)
# Find GoogleTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Include directories for the AstDyn library
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/astdyn/include)

# Add source files (you can extend this list as needed)
file(GLOB_RECURSE ASTDYN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/astdyn/src/*.cpp")
list(FILTER ASTDYN_SOURCES EXCLUDE REGEX ".*astdyn_main\\.cpp$")
list(FILTER ASTDYN_SOURCES EXCLUDE REGEX ".*astdyn_convert\\.cpp$")
add_library(astdyn STATIC ${ASTDYN_SOURCES})

target_link_libraries(astdyn Eigen3::Eigen)

# Example executable (propagation demo)
add_executable(propagate_example astdyn/examples/propagate_79148.cpp)

# Enable testing
enable_testing()

# Add test executable for AsteroidFitter
add_executable(test_asteroid_fitter astdyn/tests/test_asteroid_fitter.cpp)

target_link_libraries(test_asteroid_fitter astdyn Eigen3::Eigen GTest::GTest GTest::Main nlohmann_json::nlohmann_json)

add_test(NAME AsteroidFitterTest COMMAND test_asteroid_fitter)

# Add test executable for AsteroidFitter with Horizons comparison
add_executable(test_asteroid_fitter_horizons astdyn/tests/test_asteroid_fitter_horizons.cpp)
target_link_libraries(test_asteroid_fitter_horizons astdyn Eigen3::Eigen GTest::GTest GTest::Main nlohmann_json::nlohmann_json)
add_test(NAME AsteroidFitterHorizonsTest COMMAND test_asteroid_fitter_horizons)

# Benchmark 79450
add_executable(benchmark_79450 astdyn/tests/benchmark_79450.cpp)
target_link_libraries(benchmark_79450 astdyn Eigen3::Eigen nlohmann_json::nlohmann_json)

target_link_libraries(propagate_example astdyn Eigen3::Eigen)
