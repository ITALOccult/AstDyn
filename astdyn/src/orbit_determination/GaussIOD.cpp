/**
 * @file GaussIOD.cpp
 * @brief Implementation of Gauss IOD method
 * @author ITALOccult AstDyn Team
 * @date 2025-11-25
 */

#include "astdyn/orbit_determination/GaussIOD.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/coordinates/KeplerianElements.hpp"
#include "astdyn/time/TimeScale.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace astdyn::orbit_determination {

using namespace astdyn::constants;
using namespace astdyn::coordinates;

// ============================================================================
// GaussIODResult
// ============================================================================

void GaussIODResult::print_summary() const {
    std::cout << "\n=== Gauss IOD Result ===\n";
    std::cout << "Success: " << (success ? "Yes" : "No") << "\n";
    
    if (!success) {
        std::cout << "Error: " << error_message << "\n";
        return;
    }
    
    std::cout << "Epoch: MJD " << epoch_mjd_tdb << "\n";
    std::cout << "Observations used: #" << obs_index_1 << ", #" << obs_index_2 
              << ", #" << obs_index_3 << "\n";
    std::cout << "Iterations: " << iterations << "\n";
    std::cout << "Slant ranges: " << slant_range_1 << ", " << slant_range_2 
              << ", " << slant_range_3 << " AU\n";
    
    // Display position and velocity
    std::cout << "\nCartesian State:\n";
    std::cout << "  r = [" << state.position[0] << ", " << state.position[1] << ", " << state.position[2] << "] AU\n";
    std::cout << "  v = [" << state.velocity[0] << ", " << state.velocity[1] << ", " << state.velocity[2] << "] AU/day\n";
}

// ============================================================================
// GaussIOD
// ============================================================================

GaussIOD::GaussIOD(const GaussIODSettings& settings)
    : settings_(settings) {}

GaussIODResult GaussIOD::compute(
    const std::vector<observations::OpticalObservation>& observations) {
    
    GaussIODResult result;
    result.success = false;
    
    // Validate input
    if (observations.size() < 3) {
        result.error_message = "At least 3 observations required";
        return result;
    }
    
    // Select three optimal observations
    auto selection = select_observations(observations);
    if (!selection) {
        result.error_message = "Could not select suitable observations";
        return result;
    }
    
    auto [idx1, idx2, idx3] = *selection;
    result.obs_index_1 = idx1;
    result.obs_index_2 = idx2;
    result.obs_index_3 = idx3;
    
    if (settings_.verbose) {
        std::cout << "Selected observations: #" << idx1 << ", #" << idx2 
                  << ", #" << idx3 << "\n";
    }
    
    // Compute orbit from selected observations
    return compute_from_three(
        observations[idx1],
        observations[idx2],
        observations[idx3]);
}

GaussIODResult GaussIOD::compute_from_three(
    const observations::OpticalObservation& obs1,
    const observations::OpticalObservation& obs2,
    const observations::OpticalObservation& obs3) {
    
    GaussIODResult result;
    result.success = false;
    
    // Time intervals (in days)
    double t1 = obs1.mjd_utc;
    double t2 = obs2.mjd_utc;
    double t3 = obs3.mjd_utc;
    
    double tau1 = t1 - t2;  // Negative
    double tau3 = t3 - t2;  // Positive
    double tau = tau3 - tau1;
    
    if (settings_.verbose) {
        std::cout << "Time intervals: tau1=" << tau1 << ", tau3=" << tau3 
                  << ", tau=" << tau << " days\n";
    }
    
    // Convert UTC to TDB
    double t1_tdb = time::utc_to_tdb(t1);
    double t2_tdb = time::utc_to_tdb(t2);
    double t3_tdb = time::utc_to_tdb(t3);
    
    result.epoch_mjd_tdb = t2_tdb;
    
    // Get Earth positions at observation times (heliocentric, ecliptic J2000)
    auto earth1 = ephemeris::PlanetaryEphemeris::getState(
        ephemeris::CelestialBody::EARTH, t1_tdb);
    auto earth2 = ephemeris::PlanetaryEphemeris::getState(
        ephemeris::CelestialBody::EARTH, t2_tdb);
    auto earth3 = ephemeris::PlanetaryEphemeris::getState(
        ephemeris::CelestialBody::EARTH, t3_tdb);
    
    Vector3d R1 = earth1.position();
    Vector3d R2 = earth2.position();
    Vector3d R3 = earth3.position();
    
    // Compute line-of-sight unit vectors
    Vector3d los1 = compute_line_of_sight(obs1.ra, obs1.dec);
    Vector3d los2 = compute_line_of_sight(obs2.ra, obs2.dec);
    Vector3d los3 = compute_line_of_sight(obs3.ra, obs3.dec);
    
    // Solve for slant ranges
    double rho1, rho2, rho3;
    int iterations;
    
    bool converged = solve_slant_ranges(
        tau1, tau3, los1, los2, los3, R1, R2, R3,
        rho1, rho2, rho3, iterations);
    
    if (!converged) {
        result.error_message = "Slant range solution did not converge after " 
                             + std::to_string(iterations) + " iterations";
        return result;
    }
    
    result.slant_range_1 = rho1;
    result.slant_range_2 = rho2;
    result.slant_range_3 = rho3;
    result.iterations = iterations;
    
    // Compute heliocentric position at middle observation
    Vector3d r2 = R2 + rho2 * los2;
    
    // Estimate velocity using f,g functions
    // First get approximate r1 and r3
    Vector3d r1 = R1 + rho1 * los1;
    Vector3d r3 = R3 + rho3 * los3;
    
    // Use Gibbs method or simple finite differences for velocity
    // Simple approach: v2 ≈ (r3 - r1) / (t3 - t1)
    Vector3d v2_approx = (r3 - r1) / tau;
    
    // Refine using f,g series
    // For now use approximation - full implementation would iterate
    Vector3d v2 = v2_approx;
    
    // Store result
    result.state.position = r2;
    result.state.velocity = v2;
    result.state.epoch_mjd_tdb = t2_tdb;
    result.success = true;
    
    if (settings_.verbose) {
        std::cout << "Converged in " << iterations << " iterations\n";
        std::cout << "Position: [" << r2[0] << ", " << r2[1] << ", " << r2[2] << "] AU\n";
        std::cout << "Velocity: [" << v2[0] << ", " << v2[1] << ", " << v2[2] << "] AU/day\n";
    }
    
    return result;
}

// ============================================================================
// Private methods
// ============================================================================

std::optional<std::array<int, 3>> GaussIOD::select_observations(
    const std::vector<observations::OpticalObservation>& observations) {
    
    int n = observations.size();
    
    // Strategy: select middle observation and find two others with good separation
    int idx2 = n / 2;  // Middle observation
    
    // Find observation before middle with good separation
    int idx1 = -1;
    for (int i = idx2 - 1; i >= 0; --i) {
        double sep = observations[idx2].mjd_utc - observations[i].mjd_utc;
        if (sep >= settings_.min_separation_days && sep <= settings_.max_separation_days) {
            idx1 = i;
            break;
        }
    }
    
    // Find observation after middle with good separation
    int idx3 = -1;
    for (int i = idx2 + 1; i < n; ++i) {
        double sep = observations[i].mjd_utc - observations[idx2].mjd_utc;
        if (sep >= settings_.min_separation_days && sep <= settings_.max_separation_days) {
            idx3 = i;
            break;
        }
    }
    
    if (idx1 < 0 || idx3 < 0) {
        return std::nullopt;
    }
    
    return std::array<int, 3>{idx1, idx2, idx3};
}

Vector3d GaussIOD::compute_line_of_sight(double ra, double dec) const {
    // Convert RA/Dec to unit vector in equatorial frame
    double cos_dec = std::cos(dec);
    return Vector3d{
        cos_dec * std::cos(ra),
        cos_dec * std::sin(ra),
        std::sin(dec)
    };
}

bool GaussIOD::solve_slant_ranges(
    double tau1, double tau3,
    const Vector3d& los1, const Vector3d& los2, const Vector3d& los3,
    const Vector3d& R1, const Vector3d& R2, const Vector3d& R3,
    double& rho1, double& rho2, double& rho3,
    int& iterations) {
    
    // Gauss method with 8th order polynomial
    // Reference: Bate, Mueller & White, Ch. 5
    
    // Initial guess for rho2 (distance to object at middle observation)
    rho2 = 1.0;  // AU
    
    double tau_sq = tau3 - tau1;
    tau_sq = tau_sq * tau_sq;
    
    // Scalar triple products
    auto D0 = los1.dot(los2.cross(los3));
    if (std::abs(D0) < 1e-10) {
        iterations = 0;
        return false;  // Observations are coplanar
    }
    
    // Precompute some values
    auto D11 = R1.dot(los2.cross(los3)) / D0;
    auto D21 = R2.dot(los2.cross(los3)) / D0;
    auto D31 = R3.dot(los2.cross(los3)) / D0;
    
    auto D12 = R1.dot(los1.cross(los3)) / D0;
    auto D22 = R2.dot(los1.cross(los3)) / D0;
    auto D32 = R3.dot(los1.cross(los3)) / D0;
    
    auto D13 = R1.dot(los1.cross(los2)) / D0;
    auto D23 = R2.dot(los1.cross(los2)) / D0;
    auto D33 = R3.dot(los1.cross(los2)) / D0;
    
    // Iteration loop
    for (iterations = 0; iterations < settings_.max_iterations; ++iterations) {
        // Approximate position at middle time
        Vector3d r2 = R2 + rho2 * los2;
        double r2_mag = r2.norm();
        
        if (r2_mag < 0.01) {  // Too close to Sun
            return false;
        }
        
        // Compute f,g series coefficients (2nd order approximation)
        double u = GMS / (r2_mag * r2_mag * r2_mag);
        
        double f1 = 1.0 - 0.5 * u * tau1 * tau1;
        double g1 = tau1 - (1.0/6.0) * u * tau1 * tau1 * tau1;
        
        double f3 = 1.0 - 0.5 * u * tau3 * tau3;
        double g3 = tau3 - (1.0/6.0) * u * tau3 * tau3 * tau3;
        
        // Compute c1 and c3
        double c1 = g3 / (f1 * g3 - f3 * g1);
        double c3 = -g1 / (f1 * g3 - f3 * g1);
        
        // New rho values
        double rho1_new = (1.0/D0) * (D11 + c1 * D21 + c3 * D31);
        double rho2_new = (1.0/D0) * (D12 + c1 * D22 + c3 * D32);
        double rho3_new = (1.0/D0) * (D13 + c1 * D23 + c3 * D33);
        
        // Check for negative ranges
        if (rho1_new < 0 || rho2_new < 0 || rho3_new < 0) {
            return false;
        }
        
        // Check convergence
        double delta = std::abs(rho2_new - rho2);
        
        rho1 = rho1_new;
        rho2 = rho2_new;
        rho3 = rho3_new;
        
        if (delta < settings_.tolerance) {
            return true;  // Converged
        }
    }
    
    return false;  // Did not converge
}

std::pair<double, double> GaussIOD::compute_f_g_coefficients(
    const Vector3d& r, const Vector3d& v, double dt, double mu) const {
    
    double r_mag = r.norm();
    double v_mag = v.norm();
    double rdotv = r.dot(v);
    
    double alpha = 2.0/r_mag - v_mag*v_mag/mu;  // 1/a
    
    if (alpha > 1e-6) {  // Elliptic
        double E0 = std::acos(1.0 - r_mag * alpha);
        double n = std::sqrt(mu * alpha * alpha * alpha);
        double E = E0 + n * dt;
        
        double f = 1.0 - (mu / r_mag) * (1.0 - std::cos(E - E0));
        double g = dt - (1.0 / n) * (E - E0 - std::sin(E - E0));
        
        return {f, g};
    } else {
        // Parabolic/hyperbolic - use series expansion
        double u = mu / (r_mag * r_mag * r_mag);
        double f = 1.0 - 0.5 * u * dt * dt;
        double g = dt - (1.0/6.0) * u * dt * dt * dt;
        return {f, g};
    }
}

} // namespace astdyn::orbit_determination
