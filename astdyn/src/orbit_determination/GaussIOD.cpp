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
#include "astdyn/orbit_determination/Residuals.hpp"
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
    
    std::cout << "Epoch: MJD " << epoch.mjd() << "\n";
    std::cout << "Observations used: #" << obs_index_1 << ", #" << obs_index_2 
              << ", #" << obs_index_3 << "\n";
    std::cout << "Iterations: " << iterations << "\n";
    std::cout << "Slant ranges: " << slant_range_1 << ", " << slant_range_2 
              << ", " << slant_range_3 << " AU\n";
    
    // Display position and velocity
    std::cout << "\nCartesian State:\n";
    std::cout << "  r = [" << state.position.x_si() << ", " << state.position.y_si() << ", " << state.position.z_si() << "] m\n";
    std::cout << "  v = [" << state.velocity.x_si() << ", " << state.velocity.y_si() << ", " << state.velocity.z_si() << "] m/s\n";
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
    
    if (observations.size() < 3) {
        result.error_message = "At least 3 observations required";
        return result;
    }
    
    auto selection = select_observations(observations);
    if (!selection) {
        result.error_message = "Could not select suitable observations";
        return result;
    }
    
    auto [idx1, idx2, idx3] = *selection;
    result.obs_index_1 = idx1;
    result.obs_index_2 = idx2;
    result.obs_index_3 = idx3;
    
    return compute_from_three(observations[idx1], observations[idx2], observations[idx3]);
}

GaussIODResult GaussIOD::compute_from_three(
    const observations::OpticalObservation& obs1,
    const observations::OpticalObservation& obs2,
    const observations::OpticalObservation& obs3) {
    
    GaussIODResult result;
    result.success = false;
    
    // Time intervals (in days)
    double t1 = obs1.time.mjd();
    double t2 = obs2.time.mjd();
    double t3 = obs3.time.mjd();
    
    double tau1 = t1 - t2; 
    double tau3 = t3 - t2; 
    double tau = tau3 - tau1;
    
    // Convert to TDB
    auto t1_tdb = astdyn::time::to_tdb(obs1.time);
    auto t2_tdb = astdyn::time::to_tdb(obs2.time);
    auto t3_tdb = astdyn::time::to_tdb(obs3.time);
    
    result.epoch = t2_tdb;
    // result.state is CartesianStateTyped<core::GCRF>, it will be initialized later via from_si
    
    // Get Earth positions
    auto earth1 = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t1_tdb);
    auto earth2 = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t2_tdb);
    auto earth3 = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t3_tdb);
    
    // State::position() returns Vector3d, we need math::Vector3
    math::Vector3<core::GCRF, physics::Distance> R1 = math::Vector3<core::GCRF, physics::Distance>::from_si(earth1.position().x(), earth1.position().y(), earth1.position().z());
    math::Vector3<core::GCRF, physics::Distance> R2 = math::Vector3<core::GCRF, physics::Distance>::from_si(earth2.position().x(), earth2.position().y(), earth2.position().z());
    math::Vector3<core::GCRF, physics::Distance> R3 = math::Vector3<core::GCRF, physics::Distance>::from_si(earth3.position().x(), earth3.position().y(), earth3.position().z());
    
    // Compute line-of-sight unit vectors
    auto los1 = compute_line_of_sight(obs1.ra, obs1.dec);
    auto los2 = compute_line_of_sight(obs2.ra, obs2.dec);
    auto los3 = compute_line_of_sight(obs3.ra, obs3.dec);
    
    // Solve for slant ranges (in AU)
    double rho1, rho2, rho3;
    int iterations;
    
    bool converged = solve_slant_ranges(tau1, tau3, los1, los2, los3, R1, R2, R3, rho1, rho2, rho3, iterations);
    
    if (!converged) {
        result.error_message = "Slant range solution did not converge";
        return result;
    }
    
    result.slant_range_1 = rho1;
    result.slant_range_2 = rho2;
    result.slant_range_3 = rho3;
    result.iterations = iterations;
    
    // rho is in AU, convert to meters
    double au_to_m = constants::AU_TO_KM * 1000.0;
    
    // Compute heliocentric position at middle observation
    auto r2 = R2 + (los2 * (rho2 * au_to_m));
    
    // Estimate velocity
    auto r1 = R1 + (los1 * (rho1 * au_to_m));
    auto r3 = R3 + (los3 * (rho3 * au_to_m));
    
    auto v2 = (r3 - r1) / (tau * 86400.0);
    
    result.state = physics::CartesianStateTyped<core::GCRF>::from_si(
        t2_tdb,
        r2.x_si(), r2.y_si(), r2.z_si(),
        v2.x_si(), v2.y_si(), v2.z_si(),
        constants::GM_SUN * 1e9
    );
    result.success = true;
    
    return result;
}

std::optional<std::array<int, 3>> GaussIOD::select_observations(
    const std::vector<observations::OpticalObservation>& observations) {
    
    int n = observations.size();
    int idx2 = n / 2;
    
    int idx1 = -1;
    for (int i = idx2 - 1; i >= 0; --i) {
        double sep = observations[idx2].time.mjd() - observations[i].time.mjd();
        if (sep >= settings_.min_separation_days) {
            idx1 = i;
            break;
        }
    }
    
    int idx3 = -1;
    for (int i = idx2 + 1; i < n; ++i) {
        double sep = observations[i].time.mjd() - observations[idx2].time.mjd();
        if (sep >= settings_.min_separation_days) {
            idx3 = i;
            break;
        }
    }
    
    if (idx1 < 0 || idx3 < 0) return std::nullopt;
    return std::array<int, 3>{idx1, idx2, idx3};
}

math::Vector3<core::GCRF, physics::Distance> GaussIOD::compute_line_of_sight(double ra, double dec) const {
    double cos_dec = std::cos(dec);
    return math::Vector3<core::GCRF, physics::Distance>::from_si(
        cos_dec * std::cos(ra),
        cos_dec * std::sin(ra),
        std::sin(dec)
    );
}

bool GaussIOD::solve_slant_ranges(
    double tau1, double tau3,
    const math::Vector3<core::GCRF, physics::Distance>& los1, 
    const math::Vector3<core::GCRF, physics::Distance>& los2, 
    const math::Vector3<core::GCRF, physics::Distance>& los3,
    const math::Vector3<core::GCRF, physics::Distance>& R1, 
    const math::Vector3<core::GCRF, physics::Distance>& R2, 
    const math::Vector3<core::GCRF, physics::Distance>& R3,
    double& rho1, double& rho2, double& rho3,
    int& iterations) {
    
    // Everything here in AU for classical Gauss math
    double au_m = constants::AU_TO_KM * 1000.0;
    auto l1 = los1.to_eigen_si();
    auto l2 = los2.to_eigen_si();
    auto l3 = los3.to_eigen_si();
    auto r1_au = R1.to_eigen_si() / au_m;
    auto r2_au = R2.to_eigen_si() / au_m;
    auto r3_au = R3.to_eigen_si() / au_m;

    rho2 = 1.0; 
    
    auto D0 = l1.dot(l2.cross(l3));
    if (std::abs(D0) < 1e-10) return false;
    
    auto D11 = r1_au.dot(l2.cross(l3)) / D0;
    auto D21 = r2_au.dot(l2.cross(l3)) / D0;
    auto D31 = r3_au.dot(l2.cross(l3)) / D0;
    auto D12 = r1_au.dot(l1.cross(l3)) / D0;
    auto D22 = r2_au.dot(l1.cross(l3)) / D0;
    auto D32 = r3_au.dot(l1.cross(l3)) / D0;
    auto D13 = r1_au.dot(l1.cross(l2)) / D0;
    auto D23 = r2_au.dot(l1.cross(l2)) / D0;
    auto D33 = r3_au.dot(l1.cross(l2)) / D0;
    
    for (iterations = 0; iterations < settings_.max_iterations; ++iterations) {
        auto r2_vec = r2_au + rho2 * l2;
        double r2_mag = r2_vec.norm();
        if (r2_mag < 0.01) return false;
        
        double u = GMS / (r2_mag * r2_mag * r2_mag);
        double f1 = 1.0 - 0.5 * u * tau1 * tau1;
        double g1 = tau1 - (1.0/6.0) * u * tau1 * tau1 * tau1;
        double f3 = 1.0 - 0.5 * u * tau3 * tau3;
        double g3 = tau3 - (1.0/6.0) * u * tau3 * tau3 * tau3;
        
        double det = f1 * g3 - f3 * g1;
        if (std::abs(det) < 1e-15) return false;
        double c1 = g3 / det;
        double c3 = -g1 / det;
        
        double rho1_new = (1.0/D0) * (D11 + c1 * D21 + c3 * D31);
        double rho2_new = (1.0/D0) * (D12 + c1 * D22 + c3 * D32);
        double rho3_new = (1.0/D0) * (D13 + c1 * D23 + c3 * D33);
        
        if (rho1_new < 0 || rho2_new < 0 || rho3_new < 0) return false;
        
        double delta = std::abs(rho2_new - rho2);
        rho1 = rho1_new; rho2 = rho2_new; rho3 = rho3_new;
        
        if (delta < settings_.tolerance) return true;
    }
    
    return false;
}

std::pair<double, double> GaussIOD::compute_f_g_coefficients(
    const math::Vector3<core::GCRF, physics::Distance>& r, 
    const math::Vector3<core::GCRF, physics::Velocity>& v, 
    double dt, double mu) const {
    
    double r_mag = r.norm().to_m();
    // double v_mag = v.norm(); // unused
    // double alpha = 2.0/r_mag - v_mag*v_mag/mu; // unused
    double u = mu / (r_mag * r_mag * r_mag);
    double f = 1.0 - 0.5 * u * dt * dt;
    double g = dt - (1.0/6.0) * u * dt * dt * dt;
    return {f, g};
}

} // namespace astdyn::orbit_determination
