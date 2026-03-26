# ============================================================================
# AstDyn Tools
# ============================================================================

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

project(AstDynTools)

# Find Boost program_options
find_package(Boost REQUIRED COMPONENTS program_options)

# Occultation CLI Tool
add_executable(ioccultcalc ioccultcalc.cpp)
add_executable(astdyn_trajectory_export astdyn_trajectory_export.cpp)
add_executable(astdyn_chebyshev_demo astdyn_chebyshev_demo.cpp)
add_executable(astdyn_multi_chebyshev_demo astdyn_multi_chebyshev_demo.cpp)

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

target_link_libraries(astdyn_multi_chebyshev_demo PRIVATE 
    ${OCCULT_LIBS}
    Boost::program_options
)

# Output directory for tools
set_target_properties(ioccultcalc PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
set_target_properties(astdyn_trajectory_export PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
set_target_properties(astdyn_chebyshev_demo PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
set_target_properties(astdyn_multi_chebyshev_demo PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

# Installation
install(TARGETS ioccultcalc astdyn_trajectory_export astdyn_chebyshev_demo astdyn_multi_chebyshev_demo RUNTIME DESTINATION bin)
add_executable(haumea_star_search haumea_star_search.cpp)
target_link_libraries(haumea_star_search PRIVATE ${OCCULT_LIBS} Boost::program_options)
add_executable(haumea_prediction_test haumea_prediction_test.cpp)
target_link_libraries(haumea_prediction_test PRIVATE ${OCCULT_LIBS} Boost::program_options)
add_executable(haumea_mapper haumea_mapper.cpp)
target_link_libraries(haumea_mapper PRIVATE ${OCCULT_LIBS} Boost::program_options)
add_executable(haumea_system_integrator haumea_system_integrator.cpp)
target_link_libraries(haumea_system_integrator PRIVATE ${OCCULT_LIBS} Boost::program_options)
add_executable(haumea_system_final haumea_system_final.cpp)
target_link_libraries(haumea_system_final PRIVATE astdyn)

add_executable(haumea_baseline
    haumea_baseline.cpp
)
target_link_libraries(haumea_baseline PRIVATE astdyn)

add_executable(haumea_occult_finder haumea_occult_finder.cpp)
target_link_libraries(haumea_occult_finder PRIVATE ${OCCULT_LIBS} Boost::program_options)
add_executable(bsp_info bsp_info.cpp)
target_link_libraries(bsp_info PRIVATE ${OCCULT_LIBS} Boost::program_options)
add_executable(namaka_mapper namaka_mapper.cpp)
target_link_libraries(namaka_mapper PRIVATE ${OCCULT_LIBS} Boost::program_options)
