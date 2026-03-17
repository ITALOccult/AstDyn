cmake_minimum_required(VERSION 3.15)

# Example: Basic usage of AstDyn constants and types
add_executable(example_basic basic_usage.cpp)
target_link_libraries(example_basic PRIVATE astdyn)

# Example: Coordinate systems and reference frame transformations
# add_executable(example_coordinates coordinates_example.cpp)
# target_link_libraries(example_coordinates PRIVATE astdyn)

# Example: Planetary ephemeris calculations
# add_executable(example_ephemeris ephemeris_example.cpp)
# target_link_libraries(example_ephemeris PRIVATE astdyn)

# Example: Propagation methods benchmark
# add_executable(propagation_benchmark propagation_benchmark.cpp)
# target_link_libraries(propagation_benchmark PRIVATE astdyn)

# Example: Propagation diagnostic tool
# add_executable(propagation_diagnostic propagation_diagnostic.cpp)
# target_link_libraries(propagation_diagnostic PRIVATE astdyn)

# Example: AstDyS orbit fitter - simple workflow
# add_executable(example_astdys_fitter example_astdys_fitter.cpp)
# target_link_libraries(example_astdys_fitter PRIVATE astdyn)

# Example: Simple Pompeja orbit fitting
# add_executable(simple_pompeja_fit simple_pompeja_fit.cpp)
# target_link_libraries(simple_pompeja_fit PRIVATE astdyn)

# Fitorb API Demo
add_executable(fitorb_api_demo fitorb_api_demo.cpp)
target_link_libraries(fitorb_api_demo PRIVATE astdyn)

# Example: Configurable parsers system
# add_executable(example_configurable_parsers example_configurable_parsers.cpp)
# target_link_libraries(example_configurable_parsers PRIVATE astdyn)

# Example: Propagate 79148 for verification
# add_executable(propagate_79148 propagate_79148.cpp)
# target_link_libraries(propagate_79148 PRIVATE astdyn)

# Example: Generic propagation example for specific test case
# add_executable(propagation_example propagation_example.cpp)
# target_link_libraries(propagation_example PRIVATE astdyn)

add_executable(fitting_comparison_34713 fitting_comparison_34713.cpp)
target_link_libraries(fitting_comparison_34713 PRIVATE astdyn)

# add_executable(check_de441_positions check_de441_positions.cpp)
# target_link_libraries(check_de441_positions PRIVATE astdyn)

# JPL Horizons API Query Example
add_executable(horizons_query horizons_query.cpp)
target_link_libraries(horizons_query PRIVATE astdyn)

# Test asteroid 79518 for specific date
add_executable(test_79518 test_79518.cpp)
target_link_libraries(test_79518 PRIVATE astdyn)

# Benchmark Integrators and Models (Moved to /benchmarks)
# add_executable(benchmark_integrators benchmark_integrators.cpp)
# target_link_libraries(benchmark_integrators PRIVATE astdyn)

# Gaia Occultation Search Test
add_executable(gaia_occultation_test gaia_occultation_test.cpp)
target_link_libraries(gaia_occultation_test PRIVATE astdyn Eigen3::Eigen)

add_executable(full_od_workflow full_od_workflow.cpp)
target_link_libraries(full_od_workflow PRIVATE astdyn Eigen3::Eigen)

add_executable(critical_validation critical_validation.cpp)
target_link_libraries(critical_validation PRIVATE astdyn Eigen3::Eigen)

add_executable(debug_hzn debug_hzn.cpp)
target_link_libraries(debug_hzn PRIVATE astdyn curl)

add_executable(debug_engine debug_engine.cpp)
target_link_libraries(debug_engine PRIVATE astdyn curl Eigen3::Eigen)

add_executable(scientific_validation scientific_validation.cpp)
target_link_libraries(scientific_validation PRIVATE astdyn Eigen3::Eigen curl)

# AAS Propagation Demo
add_executable(aas_propagation_demo aas_propagation_demo.cpp)
target_link_libraries(aas_propagation_demo PRIVATE astdyn Eigen3::Eigen)

# Install examples (optional)
install(TARGETS example_basic
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/examples
)
