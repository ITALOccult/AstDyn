# ============================================================================
# AstDyn Tools
# ============================================================================

project(AstDynTools)

# Find Boost program_options
find_package(Boost REQUIRED COMPONENTS program_options)

# Occultation CLI Tool
add_executable(ioccultcalc ioccultcalc.cpp)
add_executable(astdyn_trajectory_export astdyn_trajectory_export.cpp)
add_executable(astdyn_chebyshev_demo astdyn_chebyshev_demo.cpp)

set(OCCULT_LIBS astdyn)
if(ASTDYN_BUILD_CATALOG)
    list(APPEND OCCULT_LIBS astdyn_catalog)
endif()

target_link_libraries(ioccultcalc PRIVATE 
    ${OCCULT_LIBS}
    Boost::program_options
)

target_link_libraries(astdyn_trajectory_export PRIVATE 
    ${OCCULT_LIBS}
    Boost::program_options
)

target_link_libraries(astdyn_chebyshev_demo PRIVATE 
    ${OCCULT_LIBS}
    Boost::program_options
)

# Output directory for tools
set_target_properties(ioccultcalc PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
set_target_properties(astdyn_trajectory_export PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
set_target_properties(astdyn_chebyshev_demo PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

# Installation
install(TARGETS ioccultcalc astdyn_trajectory_export astdyn_chebyshev_demo RUNTIME DESTINATION bin)
