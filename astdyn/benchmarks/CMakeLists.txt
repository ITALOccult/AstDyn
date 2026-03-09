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

# Propagation Benchmarks
add_executable(benchmark_extended src/propagation/benchmark_extended.cpp)
target_link_libraries(benchmark_extended PRIVATE astdyn)

add_executable(propagation_benchmark src/propagation/propagation_benchmark.cpp)
target_link_libraries(propagation_benchmark PRIVATE astdyn)

# Output directory for benchmarks executables
set_target_properties(
    benchmark_aas 
    benchmark_multi_integrator 
    benchmark_integrators
    od_benchmark
    benchmark_json_api
    benchmark_extended
    propagation_benchmark
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
)
