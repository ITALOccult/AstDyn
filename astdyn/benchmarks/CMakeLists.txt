# ============================================================================
# AstDyn Benchmarks
# ============================================================================

project(AstDynBenchmarks)

# Integrator Benchmarks
add_executable(benchmark_aas src/integrators/benchmark_aas.cpp)
target_link_libraries(benchmark_aas PRIVATE astdyn)

add_executable(benchmark_multi_integrator src/integrators/benchmark_multi_integrator.cpp)
target_link_libraries(benchmark_multi_integrator PRIVATE astdyn)

add_executable(benchmark_integrators src/integrators/benchmark_integrators.cpp)
target_link_libraries(benchmark_integrators PRIVATE astdyn)

# Orbit Determination Benchmarks
add_executable(od_benchmark src/orbit_determination/od_benchmark.cpp)
target_link_libraries(od_benchmark PRIVATE astdyn)

add_executable(benchmark_json_api src/orbit_determination/benchmark_json_api.cpp)
target_link_libraries(benchmark_json_api PRIVATE astdyn)

add_executable(stm_orbfit_validation src/orbit_determination/stm_orbfit_validation.cpp)
target_link_libraries(stm_orbfit_validation PRIVATE astdyn)

add_executable(od_hybrid_test src/orbit_determination/od_hybrid_test.cpp)
target_link_libraries(od_hybrid_test PRIVATE astdyn)

# Propagation Benchmarks
add_executable(benchmark_extended src/propagation/benchmark_extended.cpp)
target_link_libraries(benchmark_extended PRIVATE astdyn)

add_executable(propagation_benchmark src/propagation/propagation_benchmark.cpp)
target_link_libraries(propagation_benchmark PRIVATE astdyn)

# Occultation Benchmarks
set(OCCULT_LIBS astdyn)
if(ASTDYN_BUILD_CATALOG)
    list(APPEND OCCULT_LIBS astdyn_catalog)
endif()

add_executable(occultation_test_50936 src/occultation/occultation_test_50936.cpp)
target_link_libraries(occultation_test_50936 PRIVATE ${OCCULT_LIBS})

add_executable(star_check src/occultation/star_check.cpp)
target_link_libraries(star_check PRIVATE ${OCCULT_LIBS})

add_executable(diagnostic_50936 src/occultation/diagnostic_50936.cpp)
target_link_libraries(diagnostic_50936 PRIVATE astdyn)

add_executable(check_earth_vel src/occultation/check_earth_vel.cpp)
target_link_libraries(check_earth_vel PRIVATE astdyn)

add_executable(diagnostic_apparent src/occultation/diagnostic_apparent.cpp)
target_link_libraries(diagnostic_apparent PRIVATE astdyn)

add_executable(diagnostic_xml_tca src/occultation/diagnostic_xml_tca.cpp)
target_link_libraries(diagnostic_xml_tca PRIVATE astdyn)

add_executable(occultation_track_gen src/occultation/occultation_track_gen.cpp)
target_link_libraries(occultation_track_gen PRIVATE astdyn)

add_executable(occultation_vesta_2026 src/occultation/occultation_vesta_2026.cpp)
target_link_libraries(occultation_vesta_2026 PRIVATE ${OCCULT_LIBS})

add_executable(occultation_search_demo src/occultation/occultation_search_demo.cpp)
target_link_libraries(occultation_search_demo PRIVATE ${OCCULT_LIBS} nlohmann_json::nlohmann_json)

add_executable(occultation_catalog_benchmark src/occultation/occultation_catalog_benchmark.cpp)
target_link_libraries(occultation_catalog_benchmark PRIVATE ${OCCULT_LIBS} nlohmann_json::nlohmann_json)

add_executable(simple_catalog_benchmark src/occultation/simple_catalog_benchmark.cpp)
target_link_libraries(simple_catalog_benchmark PRIVATE ${OCCULT_LIBS} nlohmann_json::nlohmann_json)

add_executable(export_kml_bundle src/occultation/export_kml_bundle.cpp)
target_link_libraries(export_kml_bundle PRIVATE ${OCCULT_LIBS} nlohmann_json::nlohmann_json)

add_executable(occultation_system_demo src/occultation/occultation_system_demo.cpp)
target_link_libraries(occultation_system_demo PRIVATE ${OCCULT_LIBS} nlohmann_json::nlohmann_json)

add_executable(occultation_interamnia_2026 src/occultation/occultation_interamnia_2026.cpp)
target_link_libraries(occultation_interamnia_2026 PRIVATE ${OCCULT_LIBS} nlohmann_json::nlohmann_json)
add_executable(occultation_comparison_test src/occultation/occultation_comparison_test.cpp)
target_link_libraries(occultation_comparison_test PRIVATE ${OCCULT_LIBS})

add_executable(occultation_repro src/occultation/test_repro.cpp)
target_link_libraries(occultation_repro PRIVATE ${OCCULT_LIBS})

add_executable(test_xml_logic src/occultation/test_xml_logic.cpp)
target_link_libraries(test_xml_logic PRIVATE ${OCCULT_LIBS})

add_executable(occultation_final_comparison src/occultation/occultation_final_comparison.cpp)
target_link_libraries(occultation_final_comparison PRIVATE ${OCCULT_LIBS})

add_executable(test_reducer_sanity src/occultation/test_reducer_sanity.cpp)
target_link_libraries(test_reducer_sanity PRIVATE ${OCCULT_LIBS})

add_executable(test_manual_prop src/occultation/test_manual_prop.cpp)
target_link_libraries(test_manual_prop PRIVATE ${OCCULT_LIBS})

add_executable(horizons_check src/occultation/horizons_check.cpp)
target_link_libraries(horizons_check PRIVATE ${OCCULT_LIBS})

add_executable(debug_vectors src/occultation/debug_vectors.cpp)
target_link_libraries(debug_vectors PRIVATE ${OCCULT_LIBS})

# Output directory for benchmarks executables
set_target_properties(
    benchmark_aas 
    benchmark_multi_integrator 
    benchmark_integrators
    od_benchmark
    benchmark_json_api
    benchmark_extended
    propagation_benchmark
    stm_orbfit_validation
    od_hybrid_test
    occultation_test_50936
    diagnostic_50936
    star_check
    diagnostic_apparent
    diagnostic_xml_tca
    occultation_track_gen
    check_earth_vel
    occultation_vesta_2026
    occultation_search_demo
    occultation_catalog_benchmark
    simple_catalog_benchmark
    export_kml_bundle
    occultation_system_demo
    occultation_interamnia_2026
    occultation_comparison_test
    occultation_repro
    test_xml_logic
    occultation_final_comparison
    test_reducer_sanity
    test_manual_prop
    debug_vectors
    horizons_check
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
)

add_executable(horizons_check_vesta src/occultation/horizons_check_vesta.cpp)
target_link_libraries(horizons_check_vesta PRIVATE astdyn)
add_executable(horizons_diagnostic src/occultation/horizons_diagnostic.cpp)
target_link_libraries(horizons_diagnostic PRIVATE astdyn)
