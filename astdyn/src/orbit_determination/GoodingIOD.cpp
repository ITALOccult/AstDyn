/**
 * @file GoodingIOD.cpp
 * @brief Gooding's method for Initial Orbit Determination.
 */

#include "astdyn/orbit_determination/GoodingIOD.hpp"
#include "astdyn/math/LambertSolver.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/observations/ObservatoryDatabase.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::orbit_determination {

// ---------------------------------------------------------------------------
// Stumpff functions (inline, no dependency on private LambertSolver methods)
// ---------------------------------------------------------------------------
static double stumpff_C(double z) {
    if (z > 1e-6)        return (1.0 - std::cos(std::sqrt(z))) / z;
    else if (z < -1e-6)  return (std::cosh(std::sqrt(-z)) - 1.0) / (-z);
    else                 return 0.5 - z/24.0 + z*z/720.0;
}
static double stumpff_S(double z) {
    if (z > 1e-6) {
        double sz = std::sqrt(z);
        return (sz - std::sin(sz)) / (sz*sz*sz);
    } else if (z < -1e-6) {
        double sz = std::sqrt(-z);
        return (std::sinh(sz) - sz) / (sz*sz*sz);
    } else {
        return 1.0/6.0 - z/120.0 + z*z/5040.0;
    }
}

// ---------------------------------------------------------------------------
// Keplerian propagator (universal variables)
// r0 [AU], v0 [AU/day], dt_days → r_out [AU], using mu [AU³/day²]
// ---------------------------------------------------------------------------
static bool kepler_prop(const Eigen::Vector3d& r0, const Eigen::Vector3d& v0,
                        double dt_days, double mu, Eigen::Vector3d& r_out)
{
    double r0_mag = r0.norm();
    double sigma0 = r0.dot(v0) / std::sqrt(mu);
    double alpha = 2.0 / r0_mag - v0.squaredNorm() / mu; // 1/a

    double chi;
    if (alpha > 1e-9) {
        chi = std::sqrt(mu) * dt_days * alpha;
    } else if (alpha < -1e-9) {
        double a = 1.0 / alpha;
        double denom = r0.dot(v0) + std::copysign(std::sqrt(-mu * a), dt_days) * (1.0 - r0_mag * alpha);
        chi = std::copysign(1.0, dt_days) * std::sqrt(-a) * std::log(std::max(1e-300, -2.0 * mu * alpha * dt_days / denom));
    } else {
        chi = std::sqrt(mu) * dt_days / r0_mag;
    }

    for (int iter = 0; iter < 50; ++iter) {
        double psi = chi * chi * alpha;
        double c = stumpff_C(psi), s = stumpff_S(psi);
        double r_mag = sigma0*chi*(1.0 - psi*s) + (1.0 - r0_mag*alpha)*chi*chi*c + r0_mag;
        double t_pred = (sigma0*chi*chi*c + (1.0 - r0_mag*alpha)*chi*chi*chi*s + r0_mag*chi) / std::sqrt(mu);
        double dchi = (dt_days - t_pred) * std::sqrt(mu) / r_mag;
        chi += dchi;
        if (std::abs(dchi) < 1e-12) break;
    }

    double psi = chi * chi * alpha;
    double c = stumpff_C(psi), s = stumpff_S(psi);
    double r_mag = sigma0*chi*(1.0 - psi*s) + (1.0 - r0_mag*alpha)*chi*chi*c + r0_mag;
    if (r_mag < 1e-10) return false;

    double f = 1.0 - chi*chi*c / r0_mag;
    double g = dt_days - chi*chi*chi*s / std::sqrt(mu);
    r_out = f * r0 + g * v0;
    return true;
}

GoodingIOD::GoodingIOD(const Settings& settings) : settings_(settings) {}

GoodingIODResult GoodingIOD::compute(
    const observations::OpticalObservation& obs1,
    const observations::OpticalObservation& obs2,
    const observations::OpticalObservation& obs3,
    double rho1_guess_au,
    double rho3_guess_au)
{
    GoodingIODResult result;
    
    // 1. Prepare data
    time::EpochTDB t1 = time::to_tdb(obs1.time);
    time::EpochTDB t2 = time::to_tdb(obs2.time);
    time::EpochTDB t3 = time::to_tdb(obs3.time);

    auto compute_los = [](astrometry::RightAscension ra, astrometry::Declination dec) {
        double r = ra.to_rad();
        double d = dec.to_rad();
        double cra = std::cos(r);
        double sra = std::sin(r);
        double cdc = std::cos(d);
        double sdc = std::sin(d);
        // Line of sight is dimensionless unit vector
        return Eigen::Vector3d(cdc * cra, cdc * sra, sdc);
    };

    Eigen::Vector3d L1 = compute_los(obs1.ra, obs1.dec);
    Eigen::Vector3d L2 = compute_los(obs2.ra, obs2.dec);
    Eigen::Vector3d L3 = compute_los(obs3.ra, obs3.dec);

    // Get Earth positions
    auto earth1 = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t1);
    auto earth2 = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t2);
    auto earth3 = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t3);

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
        std::cout << "  [DEBUG] L1: " << L1.transpose() << std::endl;
        std::cout << "  [DEBUG] L2: " << L2.transpose() << std::endl;
        std::cout << "  [DEBUG] L3: " << L3.transpose() << std::endl;
    }

    // 3. Iteration
    double rho1_m = physics::Distance::from_au(rho1_guess_au).to_m();
    double rho3_m = physics::Distance::from_au(rho3_guess_au).to_m();
    physics::CartesianStateTyped<core::GCRF> sol_state;

    if (solve_iteration(t1, t2, t3, L1, L2, L3, R1, R2, R3, rho1_m, rho3_m, sol_state)) {
        result.success = true;
        GoodingIODResult::Solution sol;
        sol.state = sol_state;
        sol.epoch = t1;
        sol.rho1 = physics::Distance::from_si(rho1_m).to_au();
        sol.rho3 = physics::Distance::from_si(rho3_m).to_au();
        
        // Calculate rho2 from the solution state
        Eigen::Vector3d r1v = R1.to_eigen_si() + rho1_m * L1;
        Eigen::Vector3d v1v = sol_state.velocity.to_eigen_si();
        Eigen::Vector3d r2v;
        const double au_m = physics::Distance::from_au(1.0).to_m();
        kepler_prop(r1v / au_m, v1v / (au_m / 86400.0), (t2 - t1).to_days(), settings_.mu.to_au3_d2(), r2v);
        sol.rho2 = (r2v * au_m - R2.to_eigen_si()).norm() / au_m;

        sol.rms_error = 0.0; 
        result.solutions.push_back(sol);
    } else {
        result.error_message = "Gooding method failed to converge.";
    }

    return result;
}

bool GoodingIOD::solve_iteration(
    const time::EpochTDB& t1, const time::EpochTDB& t2, const time::EpochTDB& t3,
    const Eigen::Vector3d& L1, const Eigen::Vector3d& L2, const Eigen::Vector3d& L3,
    const math::Vector3<core::GCRF, physics::Distance>& R1,
    const math::Vector3<core::GCRF, physics::Distance>& R2,
    const math::Vector3<core::GCRF, physics::Distance>& R3,
    double& rho1_m, double& rho3_m,
    physics::CartesianStateTyped<core::GCRF>& final_state)
{
    const double dt13_d = (t3 - t1).to_days();
    const double dt12_d = (t2 - t1).to_days();

    // Work in AU for better numerical conditioning in the Jacobian
    const double au_to_m = physics::Distance::from_au(1.0).to_m();
    const Eigen::Vector3d R1_au = R1.to_eigen_si() / au_to_m;
    const Eigen::Vector3d R2_au = R2.to_eigen_si() / au_to_m;
    const Eigen::Vector3d R3_au = R3.to_eigen_si() / au_to_m;
    
    double rho1_au = rho1_m / au_to_m;
    double rho3_au = rho3_m / au_to_m;

    const double ra2_obs  = std::atan2(L2.y(), L2.x());
    const double dec2_obs = std::asin(std::max(-1.0, std::min(1.0, L2.z())));

    auto angular_residuals = [&](double r1_au, double r3_au,
                                  double& f_ra, double& f_dec) -> bool {
        Eigen::Vector3d r1v = R1_au + r1_au * L1;
        Eigen::Vector3d r3v = R3_au + r3_au * L3;
        
        Eigen::Vector3d v1_aupd = math::LambertSolver::solve(r1v, r3v, dt13_d, settings_.mu.to_au3_d2());
        
        Eigen::Vector3d r2p_au;
        if (!kepler_prop(r1v, v1_aupd, dt12_d, settings_.mu.to_au3_d2(), r2p_au)) return false;
        
        Eigen::Vector3d topo = r2p_au - R2_au;
        double rho2 = topo.norm();
        if (rho2 < 1e-8) return false;
        
        double ra_pred  = std::atan2(topo.y(), topo.x());
        double dec_pred = std::asin(std::max(-1.0, std::min(1.0, topo.z() / rho2)));
        
        f_ra = ra_pred - ra2_obs;
        if (f_ra >  M_PI) f_ra -= 2.0 * M_PI;
        if (f_ra < -M_PI) f_ra += 2.0 * M_PI;
        f_ra *= std::cos(dec2_obs); 
        f_dec = dec_pred - dec2_obs;
        return true;
    };

    const double eps_au = 1e-4; // 15,000 km perturbation

    bool converged = false;
    for (int outer = 0; outer < settings_.max_iterations; ++outer) {
        double f_ra, f_dec;
        if (!angular_residuals(rho1_au, rho3_au, f_ra, f_dec)) break;

        double res_norm = std::sqrt(f_ra * f_ra + f_dec * f_dec);
        if (res_norm < settings_.tolerance_rad) {
            converged = true;
            break;
        }

        // Numerical Jacobian
        double f_ra_p1, f_dec_p1, f_ra_p3, f_dec_p3;
        if (!angular_residuals(rho1_au + eps_au, rho3_au, f_ra_p1, f_dec_p1)) break;
        if (!angular_residuals(rho1_au, rho3_au + eps_au, f_ra_p3, f_dec_p3)) break;

        Eigen::Matrix2d J;
        J(0, 0) = (f_ra_p1  - f_ra)  / eps_au;
        J(1, 0) = (f_dec_p1 - f_dec) / eps_au;
        J(0, 1) = (f_ra_p3  - f_ra)  / eps_au;
        J(1, 1) = (f_dec_p3 - f_dec) / eps_au;

        double det = J.determinant();
        if (std::abs(det) < 1e-15) break;

        Eigen::Vector2d y(f_ra, f_dec);
        Eigen::Vector2d dx = -(J.colPivHouseholderQr().solve(y));
        
        // Damped step
        double max_step = 0.5 * std::max(rho1_au, rho3_au);
        double step_norm = dx.norm();
        double step_scale = (step_norm > max_step) ? (max_step / step_norm) : 1.0;
        
        rho1_au += step_scale * dx[0];
        rho3_au += step_scale * dx[1];
        
        // Safety: stay positive and within reasonable solar system bounds
        rho1_au = std::max(1e-4, std::min(rho1_au, 100.0));
        rho3_au = std::max(1e-4, std::min(rho3_au, 100.0));
    }

    if (!converged) return false;
    
    // Convert back to meters for output
    rho1_m = rho1_au * au_to_m;
    rho3_m = rho3_au * au_to_m;
              
    Eigen::Vector3d r1v_au = R1_au + rho1_au * L1;
    Eigen::Vector3d r3v_au = R3_au + rho3_au * L3;
    
    Eigen::Vector3d v1_aupd = math::LambertSolver::solve(r1v_au, r3v_au, dt13_d, settings_.mu.to_au3_d2());
    Eigen::Vector3d v1_si = v1_aupd * physics::Velocity::from_au_d(1.0).to_ms();

    final_state = physics::CartesianStateTyped<core::GCRF>::from_si(
        t1,
        r1v_au.x() * au_to_m, r1v_au.y() * au_to_m, r1v_au.z() * au_to_m,
        v1_si.x(), v1_si.y(), v1_si.z(),
        physics::GravitationalParameter::sun().to_m3_s2()
    );

    return true;
}

} // namespace astdyn::orbit_determination
