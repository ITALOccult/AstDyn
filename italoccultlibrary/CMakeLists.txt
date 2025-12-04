# ITALOccultLibrary CMake Configuration
# Libreria per propagazione orbitale con AstDyn e conversioni frame
# Validato con JPL Horizons: 0.0003 arcsec @ 3.2 AU

cmake_minimum_required(VERSION 3.15)
project(ITALOccultLibrary VERSION 1.0.0 LANGUAGES CXX)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Opzioni compilazione
option(BUILD_SHARED_LIBS "Build shared library" OFF)
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_EXAMPLES "Build examples" ON)

# Trova dipendenze
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Boost - necessario per AstDyn
set(Boost_NO_WARN_NEW_VERSIONS ON)
find_package(Boost 1.70 REQUIRED)

find_package(AstDyn REQUIRED)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${ASTDYN_INCLUDE_DIRS}
)

# Source files
set(ITALOCCULTLIB_SOURCES
    src/eq1_parser.cpp
    src/orbital_conversions.cpp
    src/astdyn_wrapper.cpp
    src/chebyshev_approximation.cpp
    src/chebyshev_rkf78_propagation.cpp
)

# Header files (per installazione)
set(ITALOCCULTLIB_HEADERS
    include/eq1_parser.h
    include/orbital_conversions.h
    include/astdyn_wrapper.h
    include/chebyshev_approximation.h
    include/chebyshev_rkf78_propagation.h
)

# Crea libreria
add_library(italoccultlib ${ITALOCCULTLIB_SOURCES})

# Alias per compatibilità
add_library(ITALOccultLibrary::italoccultlib ALIAS italoccultlib)

# Proprietà target
target_include_directories(italoccultlib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(italoccultlib PUBLIC
    Eigen3::Eigen
    AstDyn::astdyn
)

# Compiler flags
target_compile_options(italoccultlib PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -pedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

# Definizione versione
target_compile_definitions(italoccultlib PUBLIC
    ITALOCCULTLIB_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    ITALOCCULTLIB_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    ITALOCCULTLIB_VERSION_PATCH=${PROJECT_VERSION_PATCH}
)

# Installazione
install(TARGETS italoccultlib
    EXPORT ITALOccultLibraryTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(FILES ${ITALOCCULTLIB_HEADERS}
    DESTINATION include/italoccultlib
)

# Export CMake targets
install(EXPORT ITALOccultLibraryTargets
    FILE ITALOccultLibraryTargets.cmake
    NAMESPACE ITALOccultLibrary::
    DESTINATION lib/cmake/ITALOccultLibrary
)

# Config file per find_package
include(CMakePackageConfigHelpers)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/ITALOccultLibraryConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/ITALOccultLibraryConfig.cmake"
    INSTALL_DESTINATION lib/cmake/ITALOccultLibrary
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/ITALOccultLibraryConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ITALOccultLibraryConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ITALOccultLibraryConfigVersion.cmake"
    DESTINATION lib/cmake/ITALOccultLibrary
)

# Build tests
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Build examples
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Messaggi info
message(STATUS "ITALOccultLibrary Configuration:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Eigen3 version: ${EIGEN3_VERSION}")
message(STATUS "  AstDyn found: ${ASTDYN_FOUND}")
message(STATUS "  Build tests: ${BUILD_TESTS}")
message(STATUS "  Build examples: ${BUILD_EXAMPLES}")
