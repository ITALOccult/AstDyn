cmake_minimum_required(VERSION 3.15)

# Test executable for Phase 1
add_executable(astdyn_tests
    test_constants.cpp
    test_types.cpp
)

target_link_libraries(astdyn_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_tests)

# Test executable for Phase 2 - Math
add_executable(astdyn_math_tests
    test_math.cpp
)

target_link_libraries(astdyn_math_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_math_tests)

# Test executable for Phase 2 - Time
add_executable(astdyn_time_tests
    test_time.cpp
)

target_link_libraries(astdyn_time_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_time_tests)

# Test executable for Phase 2 - Utils
add_executable(astdyn_utils_tests
    test_utils.cpp
)

target_link_libraries(astdyn_utils_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_utils_tests)

# Test executable for Phase 3 - Coordinates
add_executable(astdyn_coordinates_tests
    test_coordinates.cpp
)

target_link_libraries(astdyn_coordinates_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_coordinates_tests)

# Test executable for Phase 3 - Reference Frames
add_executable(astdyn_reference_frames_tests
    test_reference_frames.cpp
)

target_link_libraries(astdyn_reference_frames_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_reference_frames_tests)

# Test executable for Phase 3 - Jacobian Matrices
add_executable(astdyn_jacobian_tests
    test_jacobians.cpp
)

target_link_libraries(astdyn_jacobian_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_jacobian_tests)

# Test executable for Phase 6 - Propagation
add_executable(astdyn_propagation_tests
    test_propagation.cpp
)

target_link_libraries(astdyn_propagation_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_propagation_tests)

# Test executable for Phase 3 - Validation Against Known Data
add_executable(astdyn_validation_tests
    test_validation.cpp
)

target_link_libraries(astdyn_validation_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_validation_tests)

# ============================================================================
# Ephemeris Tests
# ============================================================================

add_executable(astdyn_ephemeris_tests
    test_ephemeris.cpp
    test_asteroid_perturbations.cpp
)

target_link_libraries(astdyn_ephemeris_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_ephemeris_tests)

# ============================================================================
# Propagator with Asteroids Integration Tests
# ============================================================================

add_executable(astdyn_propagator_asteroids_tests
    test_propagator_asteroids.cpp
)

target_link_libraries(astdyn_propagator_asteroids_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_propagator_asteroids_tests)

# ============================================================================
# OrbFit Engine Tests
# ============================================================================

add_executable(astdyn_engine_tests
    test_astdyn_engine.cpp
)

target_link_libraries(astdyn_engine_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_engine_tests)

# ============================================================================
# Validation Tests
# ============================================================================

# Ceres validation test
add_executable(validation_ceres
    validation_ceres.cpp
)

target_link_libraries(validation_ceres
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(validation_ceres)

# OrbFit configuration tests
add_executable(astdyn_config_tests
    test_astdyn_config.cpp
    ../src/io/AstDynConfig.cpp
)

target_link_libraries(astdyn_config_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_config_tests)

# ============================================================================
# RKF78 Debug Tests
# ============================================================================

add_executable(astdyn_rkf78_debug
    test_rkf78_debug.cpp
)

target_link_libraries(astdyn_rkf78_debug
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_rkf78_debug)

# ============================================================================
# Orbit Determination Tests (Phase 7)
# ============================================================================

add_executable(astdyn_orbit_determination_tests
    test_orbit_determination.cpp
)

target_link_libraries(astdyn_orbit_determination_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_orbit_determination_tests)

# ============================================================================
# Close Approach Tests (Phase 8)
# ============================================================================

add_executable(astdyn_close_approach_tests
    test_close_approach.cpp
)

target_link_libraries(astdyn_close_approach_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_close_approach_tests)

# ============================================================================
# Mean-Osculating Conversion Tests
# ============================================================================

add_executable(astdyn_mean_osculating_tests
    test_mean_osculating.cpp
)

target_link_libraries(astdyn_mean_osculating_tests
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_mean_osculating_tests)

# ============================================================================
# Pompeja Validation Test (Real Asteroid)
# ============================================================================

add_executable(astdyn_pompeja_validation
    test_pompeja_validation.cpp
)

target_link_libraries(astdyn_pompeja_validation
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_pompeja_validation)

# ============================================================================
# Pompeja Gauss IOD Test (Initial Orbit Determination)
# ============================================================================

add_executable(astdyn_pompeja_gauss_iod
    test_pompeja_gauss_iod.cpp
)

target_link_libraries(astdyn_pompeja_gauss_iod
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_pompeja_gauss_iod)

# ============================================================================
# Pompeja Full Dynamics Test (Complete Perturbation Model)
# ============================================================================

add_executable(astdyn_pompeja_full_dynamics
    test_pompeja_full_dynamics.cpp
)

target_link_libraries(astdyn_pompeja_full_dynamics
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_pompeja_full_dynamics)

# ============================================================================
# Pompeja Differential Correction Test (Fit vs OrbFit using observations)
# ============================================================================

add_executable(astdyn_pompeja_diffcorr
    test_pompeja_diffcorr_simple.cpp
)

target_link_libraries(astdyn_pompeja_diffcorr
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_pompeja_diffcorr)

# EQ1 Parser Test
# ============================================================================

add_executable(test_eq1_parser
    test_eq1_parser.cpp
)

target_link_libraries(test_eq1_parser
    PRIVATE
        astdyn
)

# Pompeja Fit Test (using 203.eq1, 203.rwo, 203.oel)
# ============================================================================

add_executable(astdyn_pompeja_fit
    test_pompeja_fit.cpp
)

target_link_libraries(astdyn_pompeja_fit
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_pompeja_fit)

# Pompeja Comparison Test (OrbFit vs AstDyn vs Horizons)
# ============================================================================

add_executable(astdyn_pompeja_comparison
    test_pompeja_comparison.cpp
)

target_link_libraries(astdyn_pompeja_comparison
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_pompeja_comparison)

# Single Residual Debug Test
# ============================================================================

add_executable(astdyn_single_residual_debug
    test_single_residual_debug.cpp
)

target_link_libraries(astdyn_single_residual_debug
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_single_residual_debug)

# ============================================================================
# AstDyS Orbit Fitter Utility Class Test
# ============================================================================

add_executable(astdyn_astdys_fitter
    test_astdys_fitter.cpp
)

target_link_libraries(astdyn_astdys_fitter
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_astdys_fitter)

# ============================================================================
# Residual Calculation Debug Test
# ============================================================================

add_executable(astdyn_residual_debug
    test_residual_debug.cpp
)

target_link_libraries(astdyn_residual_debug
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_residual_debug)

# ============================================================================
# Single Residual Detailed Debug Test
# ============================================================================

add_executable(astdyn_residual_single
    test_residual_single.cpp
)

target_link_libraries(astdyn_residual_single
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_residual_single)

# ============================================================================
# Pompeja Differential Correction Simple Test (with AstDyS data)
# ============================================================================

add_executable(astdyn_pompeja_diffcorr_simple
    test_pompeja_diffcorr_simple.cpp
)

target_link_libraries(astdyn_pompeja_diffcorr_simple
    PRIVATE
        astdyn
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(astdyn_pompeja_diffcorr_simple)

