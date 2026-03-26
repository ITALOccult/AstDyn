cmake_minimum_required(VERSION 3.10)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

# Project name
project(ASTDYN LANGUAGES CXX)

# Require C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen (header‑only) – if not found, fall back to include path
# Find Dependencies
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(GTest REQUIRED)
find_package(CURL REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(ZLIB REQUIRED)
find_package(SQLite3 REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

# Find OpenMP (macOS helper for Homebrew)
if(APPLE)
    execute_process(COMMAND brew --prefix libomp OUTPUT_VARIABLE LIBOMP_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(LIBOMP_PREFIX)
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_CXX_LIBRARIES "-L${LIBOMP_PREFIX}/lib -lomp")
        set(OpenMP_CXX_FOUND TRUE)
        # Prevent FindOpenMP from failing if it can't find the standard ones
        set(OpenMP_C_FOUND TRUE)
        set(OpenMP_FOUND TRUE)
    endif()
endif()

if(NOT OpenMP_CXX_FOUND)
    find_package(OpenMP REQUIRED)
endif()

include_directories(${GTEST_INCLUDE_DIRS} ${CURL_INCLUDE_DIRS} ${LIBXML2_INCLUDE_DIR} ${ZLIB_INCLUDE_DIRS} ${SQLite3_INCLUDE_DIRS})

# Add external libraries
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory(astdyn/external/ioc_gaialib)

# AstDyn Library
file(GLOB_RECURSE ASTDYN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/astdyn/src/*.cpp")
list(FILTER ASTDYN_SOURCES EXCLUDE REGEX ".*astdyn_main\\.cpp$")
list(FILTER ASTDYN_SOURCES EXCLUDE REGEX ".*astdyn_convert\\.cpp$")
list(FILTER ASTDYN_SOURCES EXCLUDE REGEX ".*astdyn/src/python/.*")
add_library(astdyn STATIC ${ASTDYN_SOURCES})

target_include_directories(astdyn PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/astdyn/include
    ${CMAKE_CURRENT_SOURCE_DIR}/astdyn
    ${CMAKE_CURRENT_SOURCE_DIR}/astdyn/external/ioc_gaialib/include
)

target_link_libraries(astdyn 
    Eigen3::Eigen 
    nlohmann_json::nlohmann_json 
    ioc_gaialib 
    ${CURL_LIBRARIES} 
    ${LIBXML2_LIBRARIES} 
    ${ZLIB_LIBRARIES} 
    ${SQLite3_LIBRARIES} 
    Threads::Threads
)

if(APPLE AND LIBOMP_PREFIX)
    target_link_libraries(astdyn "-L${LIBOMP_PREFIX}/lib -lomp")
    target_compile_options(astdyn PUBLIC "-Xpreprocessor" "-fopenmp" "-I${LIBOMP_PREFIX}/include")
else()
    target_link_libraries(astdyn OpenMP::OpenMP_CXX)
endif()

# Tools
add_subdirectory(astdyn/tools)

# Enable testing
enable_testing()
