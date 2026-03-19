#!/bin/bash
BIN="./astdyn/build/tools/bin/astdyn_trajectory_export"
TF=96835.0
T0=60310.0
TOL=1e-6
ASTEROIDS_ARGS="--asteroid 1 --asteroid 1566 --asteroid 2 --asteroid 99942 --asteroid 101955 --asteroid 33179"

mkdir -p results_grkn results_gl8 results_saba4 results_ias15

echo "--------------------------------------------------------"
echo "Running GRKN for all asteroids..."
$BIN $ASTEROIDS_ARGS --t0 $T0 --tf $TF --integrator GRKN --tolerance $TOL --output "results_grkn/ast"

echo "--------------------------------------------------------"
echo "Running GL8 for all asteroids..."
$BIN $ASTEROIDS_ARGS --t0 $T0 --tf $TF --integrator GL8 --tolerance $TOL --output "results_gl8/ast"

echo "--------------------------------------------------------"
echo "Running SABA4 for all asteroids..."
$BIN $ASTEROIDS_ARGS --t0 $T0 --tf $TF --integrator SABA4 --tolerance $TOL --output "results_saba4/ast"

echo "--------------------------------------------------------"
echo "Running IAS15 for all asteroids..."
$BIN $ASTEROIDS_ARGS --t0 $T0 --tf $TF --integrator IAS15 --tolerance $TOL --output "results_ias15/ast"

echo "Benchmark run completed."
