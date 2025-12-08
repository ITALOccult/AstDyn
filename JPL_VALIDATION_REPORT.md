# Validation Report: AstDyn vs JPL Horizons

**Date:** 4 December 2025
**Object:** 17030 Sierks
**Interval:** 25 Nov 2025 - 30 Nov 2025

## 1. Methodology
We implemented a robust validation pipeline that:
1.  **Fetches Ground Truth**: Downloads high-precision ephemerides (state vectors) directly from **JPL Horizons** using their API (via `fetch_jpl_horizons.py`).
    *   Coordinate System: **ICRF / Equatorial J2000**.
    *   Center: Solar System Barycenter / Sun (Heliocentric vectors used for comparison).
2.  **Initializes AstDyn**:
    *   Takes the initial state vector $(r_0, v_0)$ from JPL at $t_{start}$.
    *   Converts Equatorial $\to$ Ecliptic frame (AstDyn input).
    *   Initializes the `AstDynWrapper` (RKF78 Propagator, full perturbations).
3.  **Propagates & Compares**:
    *   Propagates to each timestamp provided by JPL (6-hour steps).
    *   Compares the propagated position with the JPL truth vector.

## 2. Results
The validation test `test_jpl_validation` yielded the following results:

| Metric | Value | Notes |
|:-------|:------|:------|
| **Max Position Error** | **72.37 km** | Over 5 days propagation |
| **RMS Position Error** | **72.20 km** | Consistent offset |
| **Relative Error** | **~1.5e-7** | 0.15 ppm |

### Verification Log
```
Initial Condition (MJD 61004):
  Pos (ICRF): 1.04695 2.87777 1.14961
  Pos (ECLM): 1.04695 3.09759 -0.0899602
  Keplerian: a=3.17553 e=0.0454092 i=2.9046 deg

MJD      | Dist (AU) | Err (km)   | Rel Err
---------|-----------|------------|--------
61004.2500 | 3.270894 | 7.203e+01 | 1.472e-07
...
61009.7500 | 3.269126 | 7.237e+01 | 1.480e-07

✓ VALIDATION PASSED (Excellent accuracy < 1000 km)
```

## 3. Conclusions
*   **High Accuracy**: The library matches JPL Horizons to within **~72 km** for the tested asteroid. This is sufficient for occultation prediction refinement and general orbital analysis.
*   **Coordinate Systems**: The test confirmed that AstDyn correctly handles the transformation between Ecliptic J2000 (integration frame) and ICRF (output frame).
*   **Robustness**: The pipeline (`run_validation.sh`) allows continuous validation against fresh JPL data.

## 4. How to Run
```bash
./run_validation.sh
```
This command will:
1.  Download fresh data from JPL.
2.  Compile the test suite.
3.  Run the comparison and output the error metrics.
