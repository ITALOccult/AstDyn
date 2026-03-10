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
#include "astdyn/observations/ObservatoryDatabase.hpp"
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
    
    // Get Earth position at each observation time (Heliocentric)
    auto earth1 = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t1_tdb);
    auto earth2 = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t2_tdb);
    auto earth3 = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t3_tdb);

    // Apply topocentric correction: add the observatory displacement from Earth's
    // center to each observer position vector. Without this, the IOD uses the
    // geocenter (~6371 km error), which causes ~100 arcsec angular errors for
    // close-approach NEOs.
    math::Vector3<core::GCRF, physics::Distance> R1 = earth1.position;
    math::Vector3<core::GCRF, physics::Distance> R2 = earth2.position;
    math::Vector3<core::GCRF, physics::Distance> R3 = earth3.position;

    const auto& obs_db = observations::ObservatoryDatabase::getInstance();
    if (auto info1 = obs_db.getObservatory(obs1.observatory_code)) {
        auto topo = info1->getPositionGCRF(obs1.time);
        R1 = R1 + math::Vector3<core::GCRF, physics::Distance>::from_si(topo.x_si(), topo.y_si(), topo.z_si());
    }
    if (auto info2 = obs_db.getObservatory(obs2.observatory_code)) {
        auto topo = info2->getPositionGCRF(obs2.time);
        R2 = R2 + math::Vector3<core::GCRF, physics::Distance>::from_si(topo.x_si(), topo.y_si(), topo.z_si());
    }
    if (auto info3 = obs_db.getObservatory(obs3.observatory_code)) {
        auto topo = info3->getPositionGCRF(obs3.time);
        R3 = R3 + math::Vector3<core::GCRF, physics::Distance>::from_si(topo.x_si(), topo.y_si(), topo.z_si());
    }
    
    if (settings_.verbose) {
        std::cout << "  [DEBUG] R1: " << R1.x_si()/1.495978707e11 << ", " << R1.y_si()/1.495978707e11 << ", " << R1.z_si()/1.495978707e11 << " AU" << std::endl;
        std::cout << "  [DEBUG] R2: " << R2.x_si()/1.495978707e11 << ", " << R2.y_si()/1.495978707e11 << ", " << R2.z_si()/1.495978707e11 << " AU" << std::endl;
        std::cout << "  [DEBUG] R3: " << R3.x_si()/1.495978707e11 << ", " << R3.y_si()/1.495978707e11 << ", " << R3.z_si()/1.495978707e11 << " AU" << std::endl;
    }

    // Compute line-of-sight unit vectors (dimensionless)
    Eigen::Vector3d los1 = compute_line_of_sight(obs1.ra, obs1.dec);
    Eigen::Vector3d los2 = compute_line_of_sight(obs2.ra, obs2.dec);
    Eigen::Vector3d los3 = compute_line_of_sight(obs3.ra, obs3.dec);
    
    if (settings_.verbose) {
        std::cout << "  [DEBUG] L1: " << los1.transpose() << std::endl;
        std::cout << "  [DEBUG] L2: " << los2.transpose() << std::endl;
        std::cout << "  [DEBUG] L3: " << los3.transpose() << std::endl;
    }
    
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

    // Use physics::Distance for radius vectors
    double d1 = physics::Distance::from_au(rho1).to_m();
    double d2 = physics::Distance::from_au(rho2).to_m();
    double d3 = physics::Distance::from_au(rho3).to_m();
    
    auto radius1 = R1 + math::Vector3<core::GCRF, physics::Distance>::from_si(los1.x() * d1, los1.y() * d1, los1.z() * d1);
    auto radius2 = R2 + math::Vector3<core::GCRF, physics::Distance>::from_si(los2.x() * d2, los2.y() * d2, los2.z() * d2);
    auto radius3 = R3 + math::Vector3<core::GCRF, physics::Distance>::from_si(los3.x() * d3, los3.y() * d3, los3.z() * d3);
    
    // Estimate velocity: v2 = (f1*r3 - f3*r1) / (g1*f3 - g3*f1) [m/s]
    // Use the values from the last iteration for best accuracy
    double r2_mag_au = radius2.norm().to_au();
    double mu_au = settings_.mu.to_au3_d2();
    double u = mu_au / (r2_mag_au * r2_mag_au * r2_mag_au);
    
    double f1 = 1.0 - 0.5 * u * tau1 * tau1;
    double g1 = tau1 - (1.0/6.0) * u * tau1 * tau1 * tau1;
    double f3 = 1.0 - 0.5 * u * tau3 * tau3;
    double g3 = tau3 - (1.0/6.0) * u * tau3 * tau3 * tau3;

    double det = f1 * g3 - f3 * g1;
    // Current radius1, radius3 are math::Vector3<GCRF, Distance> (m)
    // Velocities must be in m per SECOND
    auto v2_vec = (radius3 * f1 - radius1 * f3) / (det * 86400.0);
    
    result.state = physics::CartesianStateTyped<core::GCRF>::from_si(
        t2_tdb,
        radius2.x_si(), radius2.y_si(), radius2.z_si(),
        v2_vec.x_si(), v2_vec.y_si(), v2_vec.z_si(),
        physics::GravitationalParameter::sun().to_m3_s2()
    );
    result.success = true;
    
    return result;
}

std::optional<std::array<int, 3>> GaussIOD::select_observations(
    const std::vector<observations::OpticalObservation>& observations) {
    
    int n = observations.size();
    if (n < 3) return std::nullopt;

    // We want three observations that satisfy:
    // 1. Separation > min_separation_days
    // 2. Total arc < max_separation_days
    // 3. Middle observation should be roughly in the center of the triplet
    
    // Strategy: Search for a triplet that fits the constraints, starting from the beginning
    for (int i = 0; i < n - 2; ++i) {
        for (int k = i + 2; k < n; ++k) {
            double total_arc = observations[k].time.mjd() - observations[i].time.mjd();
            
            // For Gauss IOD, total arc > 30-60 days usually fails due to Taylor series divergence.
            // Optimal arc is 5-20 days for most LEOS/NEOs.
            if (total_arc >= 2.0 * settings_.min_separation_days && total_arc <= settings_.max_separation_days) {
                // Find a middle observation close to center
                double target_mjd = observations[i].time.mjd() + 0.5 * total_arc;
                int j_best = -1;
                double min_diff = 1e9;
                
                for (int j = i + 1; j < k; ++j) {
                    double sep1 = observations[j].time.mjd() - observations[i].time.mjd();
                    double sep2 = observations[k].time.mjd() - observations[j].time.mjd();
                    
                    if (sep1 >= settings_.min_separation_days && sep2 >= settings_.min_separation_days) {
                        double diff = std::abs(observations[j].time.mjd() - target_mjd);
                        if (diff < min_diff) {
                            min_diff = diff;
                            j_best = j;
                        }
                    }
                }
                
                if (j_best != -1) {
                    return std::array<int, 3>{i, j_best, k};
                }
            }
        }
    }
    
    return std::nullopt;
}

Eigen::Vector3d GaussIOD::compute_line_of_sight(astrometry::RightAscension ra, astrometry::Declination dec) const {
    double r = ra.to_rad();
    double d = dec.to_rad();
    return Eigen::Vector3d(
        std::cos(d) * std::cos(r),
        std::cos(d) * std::sin(r),
        std::sin(d)
    );
}

bool GaussIOD::solve_slant_ranges(
    double tau1, double tau3,
    const Eigen::Vector3d& l1, 
    const Eigen::Vector3d& l2, 
    const Eigen::Vector3d& l3,
    const math::Vector3<core::GCRF, physics::Distance>& R1, 
    const math::Vector3<core::GCRF, physics::Distance>& R2, 
    const math::Vector3<core::GCRF, physics::Distance>& R3,
    double& rho1, double& rho2, double& rho3,
    int& iterations) {
    
    // Everything here in AU for classical Gauss math
    double au_m = physics::Distance::from_au(1.0).to_m();
    auto r1_au = R1.to_eigen_si() / au_m;
    auto r2_au = R2.to_eigen_si() / au_m;
    auto r3_au = R3.to_eigen_si() / au_m;

    rho2 = 1.0; 
    
    auto D0 = l1.dot(l2.cross(l3));
    if (std::abs(D0) < 1e-8) {
        return false;
    }
    
    // D-matrix elements (raw, not divided by D0)
    auto D11 = r1_au.dot(l2.cross(l3));
    auto D21 = r2_au.dot(l2.cross(l3));
    auto D31 = r3_au.dot(l2.cross(l3));
    auto D12 = r1_au.dot(l1.cross(l3));
    auto D22 = r2_au.dot(l1.cross(l3));
    auto D32 = r3_au.dot(l1.cross(l3));
    auto D13 = r1_au.dot(l1.cross(l2));
    auto D23 = r2_au.dot(l1.cross(l2));
    auto D33 = r3_au.dot(l1.cross(l2));

    for (iterations = 0; iterations < settings_.max_iterations; ++iterations) {
        auto r2_vec = r2_au + rho2 * l2;
        double r2_mag = r2_vec.norm();
        if (r2_mag < 0.01) return false;

        double mu_au = settings_.mu.to_au3_d2();
        double u = mu_au / (r2_mag * r2_mag * r2_mag);
        double f1 = 1.0 - 0.5 * u * tau1 * tau1;
        double g1 = tau1 - (1.0/6.0) * u * tau1 * tau1 * tau1;
        double f3 = 1.0 - 0.5 * u * tau3 * tau3;
        double g3 = tau3 - (1.0/6.0) * u * tau3 * tau3 * tau3;

        double det = f1 * g3 - f3 * g1;
        if (std::abs(det) < 1e-15) return false;
        double c1 = g3 / det;
        double c3 = -g1 / det;

        // Cramer's rule: [l1 | -l2 | l3] * [c1*rho1; rho2; c3*rho3]^T = R2 - c1*R1 - c3*R3
        // det([l1|-l2|l3]) = -D0
        // rho1: c1*rho1 = det([b|-l2|l3]) / (-D0) = (D21 - c1*D11 - c3*D31) / D0
        // rho2:    rho2 = det([l1|b|l3])  / (-D0) = (D22 - c1*D12 - c3*D32) / D0
        // rho3: c3*rho3 = det([l1|-l2|b]) / (-D0) = (D23 - c1*D13 - c3*D33) / D0
        if (std::abs(c1) < 1e-15 || std::abs(c3) < 1e-15) return false;

        double rho1_new = (D21 - c1 * D11 - c3 * D31) / (c1 * D0);
        double rho2_new = (c1 * D12 - D22 + c3 * D32) / D0;
        double rho3_new = (D23 - c1 * D13 - c3 * D33) / (c3 * D0);

        double delta = std::abs(rho2_new - rho2);
        rho1 = rho1_new; rho2 = rho2_new; rho3 = rho3_new;

        if (delta < settings_.tolerance.to_au()) {
            // Only reject unphysical solutions after convergence
            if (rho1 < 0 || rho2 < 0 || rho3 < 0) return false;
            return true;
        }
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
