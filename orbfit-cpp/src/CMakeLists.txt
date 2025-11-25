# Collect all source files
file(GLOB_RECURSE ORBFIT_SOURCES
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

# Add OrbFitEngine to the library
list(APPEND ORBFIT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/OrbFitEngine.cpp)

# Create main library
add_library(orbfit ${ORBFIT_SOURCES})

# Set library properties
set_target_properties(orbfit PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/../include/orbfit/OrbFit.hpp"
)

# Link dependencies
target_link_libraries(orbfit
    PUBLIC
        Eigen3::Eigen
        Boost::filesystem
        Boost::program_options
        Boost::date_time
)

# Optionally link CSPICE
if(ORBFIT_USE_SPICE AND CSPICE_FOUND)
    target_link_libraries(orbfit PUBLIC CSPICE::CSPICE)
endif()

# Include directories
target_include_directories(orbfit
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Install library
install(TARGETS orbfit
    EXPORT OrbFitTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/orbfit
)

# Install export targets
install(EXPORT OrbFitTargets
    FILE OrbFitTargets.cmake
    NAMESPACE OrbFit::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/OrbFit
)

# ============================================================================
# Main Programs
# ============================================================================

# orbfit - Main orbit determination program
add_executable(orbfit_program
    orbfit_main.cpp
)

target_link_libraries(orbfit_program
    PRIVATE
        orbfit
)

set_target_properties(orbfit_program PROPERTIES
    OUTPUT_NAME orbfit
)

# orbfit_convert - Configuration converter tool
add_executable(orbfit_convert
    orbfit_convert.cpp
    io/OrbFitConfig.cpp
)

target_link_libraries(orbfit_convert
    PRIVATE
        orbfit
)

# Install programs
install(TARGETS orbfit_program orbfit_convert
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
