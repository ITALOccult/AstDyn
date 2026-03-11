#pragma once

/**
 * OrbFit C++ — Orbital fitting pipeline
 * ======================================
 * Implementation of the Milani & Gronchi OrbFit algorithm in C++17.
 *
 * Pipeline:
 *  1. Initial Orbit Determination (Gauss method)            → iod.h
 *  2. N-body propagation (Dormand-Prince RK45 + STM)        → integrator.h
 *  3. Differential Corrections (Gauss-Newton / LM)          → differential_corrections.h
 *  4. Uncertainty propagation & risk assessment             → uncertainty.h
 *
 * References:
 *  - Milani & Gronchi (2010) "Theory of Orbit Determination", Cambridge UP
 *  - Gauss (1809) "Theoria Motus Corporum Coelestium"
 *  - Standish (1992) "Orientation of the JPL Ephemerides", A&A 114
 *  - Dormand & Prince (1980) JCAM 6, 19–26
 */

#include "constants.h"
#include "linalg.h"
#include "orbital_elements.h"
#include "integrator.h"
#include "observations.h"
#include "iod.h"
#include "differential_corrections.h"
#include "uncertainty.h"

namespace orbfit {

// ─────────────────────────────────────────────────────────────
//  Full pipeline: IOD + DC
// ─────────────────────────────────────────────────────────────
struct PipelineResult {
    GaussSolution   iod;
    FitResult       fit;
    ConfidenceEllipsoid ellipsoid;
};

inline PipelineResult run_pipeline(const std::vector<Observation>& obs,
                                    const DCOptions& dco = {})
{
    if(obs.size() < 3)
        throw std::runtime_error("Pipeline requires at least 3 observations");

    // Step 1: Initial orbit determination
    auto [i1,i2,i3] = select_iod_obs(obs);
    GaussSolution iod = gauss_iod(obs[i1], obs[i2], obs[i3]);

    // Step 2: Differential corrections
    FitResult fit = differential_corrections(iod.elements, obs, dco);

    // Step 3: Confidence ellipsoid
    ConfidenceEllipsoid ell;
    ell.covariance  = fit.covariance;
    ell.sigma_level = 3.0;

    return {iod, fit, ell};
}

} // namespace orbfit
