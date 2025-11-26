cmake_minimum_required(VERSION 3.15)

# Example: Basic usage of AstDyn constants and types
add_executable(example_basic basic_usage.cpp)
target_link_libraries(example_basic PRIVATE astdyn)

# Example: Coordinate systems and reference frame transformations
add_executable(example_coordinates coordinates_example.cpp)
target_link_libraries(example_coordinates PRIVATE astdyn)

# Example: Planetary ephemeris calculations
add_executable(example_ephemeris ephemeris_example.cpp)
target_link_libraries(example_ephemeris PRIVATE astdyn)

# Example: Propagation methods benchmark
add_executable(propagation_benchmark propagation_benchmark.cpp)
target_link_libraries(propagation_benchmark PRIVATE astdyn)

# Example: Propagation diagnostic tool
add_executable(propagation_diagnostic propagation_diagnostic.cpp)
target_link_libraries(propagation_diagnostic PRIVATE astdyn)

# Example: AstDyS orbit fitter - simple workflow
add_executable(example_astdys_fitter example_astdys_fitter.cpp)
target_link_libraries(example_astdys_fitter PRIVATE astdyn)

# Example: Simple Pompeja orbit fitting (has API errors - disabled for now)
# add_executable(simple_pompeja_fit simple_pompeja_fit.cpp)
# target_link_libraries(simple_pompeja_fit PRIVATE astdyn)

# Example: Configurable parsers system
add_executable(example_configurable_parsers example_configurable_parsers.cpp)
target_link_libraries(example_configurable_parsers PRIVATE astdyn)

# Install examples (optional)
install(TARGETS example_basic example_coordinates example_ephemeris propagation_benchmark example_astdys_fitter example_configurable_parsers
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/examples
)
