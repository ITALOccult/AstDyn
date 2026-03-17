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

GoodingIOD::GoodingIOD(std::shared_ptr<ephemeris::PlanetaryEphemeris> ephem, const Settings& settings) 
    : settings_(settings), ephem_(ephem) {}

GoodingIODResult GoodingIOD::compute(
    const observations::OpticalObservation& obs1,
    const observations::OpticalObservation& obs2,
    const observations::OpticalObservation& obs3,
    physics::Distance rho1_guess,
    physics::Distance rho3_guess)
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
    auto earth1 = ephem_->getState(ephemeris::CelestialBody::EARTH, t1);
    auto earth2 = ephem_->getState(ephemeris::CelestialBody::EARTH, t2);
    auto earth3 = ephem_->getState(ephemeris::CelestialBody::EARTH, t3);

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
        std::cout << "  [DEBUG] R1: " << R1.x_si()/1.495978707e11 << ", " << R1.y_si()/1.495978707e11 << ", " << R1.z_si()/1.495978707e11 << " AU" << "\n";
        std::cout << "  [DEBUG] R2: " << R2.x_si()/1.495978707e11 << ", " << R2.y_si()/1.495978707e11 << ", " << R2.z_si()/1.495978707e11 << " AU" << "\n";
        std::cout << "  [DEBUG] R3: " << R3.x_si()/1.495978707e11 << ", " << R3.y_si()/1.495978707e11 << ", " << R3.z_si()/1.495978707e11 << " AU" << "\n";
        std::cout << "  [DEBUG] L1: " << L1.transpose() << "\n";
        std::cout << "  [DEBUG] L2: " << L2.transpose() << "\n";
        std::cout << "  [DEBUG] L3: " << L3.transpose() << "\n";
    }

    // 3. Iteration
    physics::Distance rho1 = rho1_guess;
    physics::Distance rho3 = rho3_guess;
    physics::CartesianStateTyped<core::GCRF> sol_state;

    if (solve_iteration(t1, t2, t3, L1, L2, L3, R1, R2, R3, rho1, rho3, sol_state)) {
        result.success = true;
        GoodingIODResult::Solution sol;
        sol.state = sol_state;
        sol.epoch = t1;
        sol.rho1 = rho1;
        sol.rho3 = rho3;
        
        // Calculate rho2 from the solution state
        Eigen::Vector3d r1v = R1.to_eigen_si() + rho1.to_m() * L1;
        Eigen::Vector3d v1v = sol_state.velocity.to_eigen_si();
        Eigen::Vector3d r2v;
        kepler_prop(r1v / constants::AU, v1v / (constants::AU / constants::SECONDS_PER_DAY), (t2 - t1).to_days(), settings_.mu.to_au3_d2(), r2v);
        sol.rho2 = physics::Distance::from_si((r2v * constants::AU - R2.to_eigen_si()).norm());

        sol.rms_error = astrometry::Angle::from_rad(0.0); 
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
    physics::Distance& rho1, physics::Distance& rho3,
    physics::CartesianStateTyped<core::GCRF>& final_state)
{
    const time::TimeDuration dt13 = (t3 - t1);
    const time::TimeDuration dt12 = (t2 - t1);
    const double dt12_d = dt12.to_days();

    // Use AU for numerical stability in the Newton solver
    const Eigen::Vector3d R1_au = R1.to_eigen_si() / constants::AU;
    const Eigen::Vector3d R2_au = R2.to_eigen_si() / constants::AU;
    const Eigen::Vector3d R3_au = R3.to_eigen_si() / constants::AU;
    
    double rho1_au = rho1.to_au();
    double rho3_au = rho3.to_au();

    const double ra2_obs  = std::atan2(L2.y(), L2.x());
    const double dec2_obs = std::asin(std::max(-1.0, std::min(1.0, L2.z())));

    auto solve_angular_residuals = [&](double current_rho1_au, double current_rho3_au,
                                  double& f_ra, double& f_dec) -> bool {
        auto rho1_dist = physics::Distance::from_au(current_rho1_au);
        auto rho3_dist = physics::Distance::from_au(current_rho3_au);

        auto r1v_full = R1 + math::Vector3<core::GCRF, physics::Distance>::from_si(L1.x() * rho1_dist.to_m(), L1.y() * rho1_dist.to_m(), L1.z() * rho1_dist.to_m());
        auto r3v_full = R3 + math::Vector3<core::GCRF, physics::Distance>::from_si(L3.x() * rho3_dist.to_m(), L3.y() * rho3_dist.to_m(), L3.z() * rho3_dist.to_m());
        
        auto v1_res = math::LambertSolver::solve(r1v_full, r3v_full, dt13, settings_.mu);
        
        Eigen::Vector3d r2p_au;
        Eigen::Vector3d r1v_au = r1v_full.to_eigen_si() / constants::AU;
        Eigen::Vector3d v1v_aupd = v1_res.to_eigen_si() / (constants::AU / constants::SECONDS_PER_DAY);

        if (!kepler_prop(r1v_au, v1v_aupd, dt12_d, settings_.mu.to_au3_d2(), r2p_au)) return false;
        
        Eigen::Vector3d topo = r2p_au - R2_au;
        double rho2 = topo.norm();
        if (rho2 < 1e-8) return false;
        
        double ra_pred  = std::atan2(topo.y(), topo.x());
        double dec_pred = std::asin(std::max(-1.0, std::min(1.0, topo.z() / rho2)));
        
        f_ra = ra_pred - ra2_obs;
        while (f_ra >  constants::PI) f_ra -= 2.0 * constants::PI;
        while (f_ra < -constants::PI) f_ra += 2.0 * constants::PI;
        f_ra *= std::cos(dec2_obs); 
        f_dec = dec_pred - dec2_obs;
        return true;
    };

    const double perturbation_au = 1e-4; 
    bool converged = false;

    for (int iter = 0; iter < settings_.max_iterations; ++iter) {
        double f_ra, f_dec;
        if (!solve_angular_residuals(rho1_au, rho3_au, f_ra, f_dec)) break;

        if (std::sqrt(f_ra * f_ra + f_dec * f_dec) < settings_.tolerance_rad) {
            converged = true;
            break;
        }

        double f_ra1, f_dec1, f_ra3, f_dec3;
        if (!solve_angular_residuals(rho1_au + perturbation_au, rho3_au, f_ra1, f_dec1)) break;
        if (!solve_angular_residuals(rho1_au, rho3_au + perturbation_au, f_ra3, f_dec3)) break;

        Eigen::Matrix2d jacobian;
        jacobian(0, 0) = (f_ra1  - f_ra)  / perturbation_au;
        jacobian(1, 0) = (f_dec1 - f_dec) / perturbation_au;
        jacobian(0, 1) = (f_ra3  - f_ra)  / perturbation_au;
        jacobian(1, 1) = (f_dec3 - f_dec) / perturbation_au;

        if (std::abs(jacobian.determinant()) < 1e-15) break;

        Eigen::Vector2d residuals(f_ra, f_dec);
        Eigen::Vector2d step = -(jacobian.colPivHouseholderQr().solve(residuals));
        
        double max_allowed = 0.5 * std::max(rho1_au, rho3_au);
        double step_magnitude = step.norm();
        double scale = (step_magnitude > max_allowed) ? (max_allowed / step_magnitude) : 1.0;
        
        rho1_au += scale * step[0];
        rho3_au += scale * step[1];
        
        rho1_au = std::max(1e-4, std::min(rho1_au, 100.0));
        rho3_au = std::max(1e-4, std::min(rho3_au, 100.0));
    }

    if (!converged) return false;
    
    Eigen::Vector3d r1v_au = R1_au + rho1_au * L1;
    Eigen::Vector3d r3v_au = R3_au + rho3_au * L3;
    
    auto r1v_typed = math::Vector3<core::GCRF, physics::Distance>::from_si(
        physics::Distance::from_au(r1v_au.x()).to_m(),
        physics::Distance::from_au(r1v_au.y()).to_m(),
        physics::Distance::from_au(r1v_au.z()).to_m());
    auto r3v_typed = math::Vector3<core::GCRF, physics::Distance>::from_si(
        physics::Distance::from_au(r3v_au.x()).to_m(),
        physics::Distance::from_au(r3v_au.y()).to_m(),
        physics::Distance::from_au(r3v_au.z()).to_m());
    auto v1_typed = math::LambertSolver::solve(r1v_typed, r3v_typed, dt13, settings_.mu);

    final_state = physics::CartesianStateTyped<core::GCRF>::from_si(
        t1, r1v_au.x() * constants::AU, r1v_au.y() * constants::AU, r1v_au.z() * constants::AU,
        v1_typed.x_si(), v1_typed.y_si(), v1_typed.z_si(),
        settings_.mu.to_m3_s2()
    );

    return true;
}

} // namespace astdyn::orbit_determination
