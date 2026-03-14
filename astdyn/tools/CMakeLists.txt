# ============================================================================
# AstDyn Tools
# ============================================================================

project(AstDynTools)

# Find Boost program_options
find_package(Boost REQUIRED COMPONENTS program_options)

# Occultation CLI Tool
add_executable(ioccultcalc ioccultcalc.cpp)

set(OCCULT_LIBS astdyn)
if(ASTDYN_BUILD_CATALOG)
    list(APPEND OCCULT_LIBS astdyn_catalog)
endif()

target_link_libraries(ioccultcalc PRIVATE 
    ${OCCULT_LIBS}
    Boost::program_options
)

# Output directory for tools
set_target_properties(ioccultcalc PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

# Installation
install(TARGETS ioccultcalc RUNTIME DESTINATION bin)
