cmake_minimum_required(VERSION 3.15)

# Unified AstDyn Unit Tests
add_executable(astdyn_unit_tests
    test_constants.cpp
    test_types.cpp
    test_time.cpp
    test_utils.cpp
    test_coordinates.cpp
    test_reference_frames.cpp
    test_jacobians.cpp
    test_propagation.cpp
    test_astdyn_engine.cpp
    test_rigorous_validation.cpp
    test_integrator_comparative.cpp
    test_occultation.cpp
    test_asteroid_perturbations.cpp
)

target_link_libraries(astdyn_unit_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)
if(ASTDYN_BUILD_CATALOG)
    target_link_libraries(astdyn_unit_tests PRIVATE astdyn_catalog)
endif()

gtest_discover_tests(astdyn_unit_tests)

# Benchmarks moved to /benchmarks directory
# add_executable(astdyn_benchmark_aas ...)
# add_executable(benchmark_multi_integrator ...)
# add_executable(benchmark_extended ...)

# Fetch Horizons Vectors Tool (Kept in tools, built here for convenience)
add_executable(fetch_horizons_vectors
    ../tools/fetch_horizons_vectors.cpp
)

target_link_libraries(fetch_horizons_vectors
    PRIVATE
        astdyn
        curl
)
