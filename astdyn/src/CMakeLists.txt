# Collect all source files
file(GLOB_RECURSE ASTDYN_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/core/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/math/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/time/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/coordinates/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/orbit/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/io/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ephemeris/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/observations/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/propagation/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/orbit_determination/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/close_approach/*.cpp
)

# Add AstDynEngine to the library
list(APPEND ASTDYN_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/AstDynEngine.cpp)

# Create main library
add_library(astdyn ${ASTDYN_SOURCES})

# Set library properties
set_target_properties(astdyn PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/../include/astdyn/AstDyn.hpp"
)

# Link dependencies
target_link_libraries(astdyn
    PUBLIC
        Eigen3::Eigen
        Boost::filesystem
        Boost::program_options
        Boost::date_time
)

# Optionally link CSPICE
if(ASTDYN_USE_SPICE AND CSPICE_FOUND)
    target_link_libraries(astdyn PUBLIC CSPICE::CSPICE)
endif()

# Include directories
target_include_directories(astdyn
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Install library
install(TARGETS astdyn
    EXPORT AstDynTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/astdyn
)

# Install export targets
install(EXPORT AstDynTargets
    FILE AstDynTargets.cmake
    NAMESPACE AstDyn::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AstDyn
)

# ============================================================================
# Main Programs
# ============================================================================

# astdyn - Main orbit determination program
add_executable(astdyn_program
    astdyn_main.cpp
)

target_link_libraries(astdyn_program
    PRIVATE
        astdyn
)

set_target_properties(astdyn_program PROPERTIES
    OUTPUT_NAME astdyn
)

# astdyn_convert - Configuration converter tool
add_executable(astdyn_convert
    astdyn_convert.cpp
    io/AstDynConfig.cpp
)

target_link_libraries(astdyn_convert
    PRIVATE
        astdyn
)

# Install programs
install(TARGETS astdyn_program astdyn_convert
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
