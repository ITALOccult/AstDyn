#!/bin/bash
# Final Integrator Benchmark Comparison Script
BIN="./astdyn/build/tools/bin/astdyn_trajectory_export"
TF=70000.0  # Approx 26 years export
T0=60310.0
TOL=1e-8
ASTEROIDS_ARGS="--asteroid 1 --asteroid 1566 --asteroid 2 --asteroid 99942 --asteroid 101955 --asteroid 33179"

INTEGRATORS=("AAS" "RKF78" "GRKN" "GL8" "SABA4" "IAS15")

for INT in "${INTEGRATORS[@]}"; do
    DIR="results_$(echo $INT | tr '[:upper:]' '[:lower:]')"
    mkdir -p $DIR
    echo "--------------------------------------------------------"
    echo "Running Benchmark for integrator: $INT"
    $BIN $ASTEROIDS_ARGS --t0 $T0 --tf $TF --integrator $INT --tolerance $TOL --output "$DIR/ast"
done

echo "Full comparison run completed."
