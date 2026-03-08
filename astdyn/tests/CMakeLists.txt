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
)

target_link_libraries(astdyn_unit_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_unit_tests)

# Benchmark AAS tool (Keep separate as it's a CLI tool, not a unit test)
add_executable(astdyn_benchmark_aas
    ../tools/benchmark_aas.cpp
)

target_link_libraries(astdyn_benchmark_aas
    PRIVATE
        astdyn
)
