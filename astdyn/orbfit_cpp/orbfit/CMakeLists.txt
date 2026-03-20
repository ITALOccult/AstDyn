cmake_minimum_required(VERSION 3.14)
project(OrbFitCpp VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -fsanitize=address")

include_directories(include)

# Demo executable
add_executable(orbfit_demo examples/demo.cpp)
target_include_directories(orbfit_demo PRIVATE include)

# Install headers
install(DIRECTORY include/ DESTINATION include/orbfit)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
