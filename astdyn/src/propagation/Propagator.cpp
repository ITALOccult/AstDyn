#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/propagation/ForceField.hpp"
#include "astdyn/propagation/AASIntegrator.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::propagation {

Propagator::Propagator(std::shared_ptr<Integrator> integrator,
                      std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris,
                      const PropagatorSettings& settings)
    : integrator_(std::move(integrator)), ephemeris_(std::move(ephemeris)), settings_(settings) {
    mat_ecl_ = coordinates::ReferenceFrame::j2000_to_ecliptic();
    setup_aas_parameters();
}

void Propagator::setup_aas_parameters() {
    auto aas = std::dynamic_pointer_cast<AASIntegrator>(integrator_);
    if (!aas) return;
    bool sun = std::abs(settings_.central_body_gm - constants::GMS) < 1e-5;
    double j2 = sun ? (settings_.include_sun_j2 ? constants::SUN_J2 : 0.0) 
                    : (settings_.include_earth_j2 ? constants::EARTH_J2 : 0.0);
    double r_eq = sun ? constants::R_SUN_AU : (constants::R_EARTH / constants::AU);
    aas->set_central_body(settings_.central_body_gm, j2, r_eq);
}

Eigen::VectorXd Propagator::compute_derivatives(time::EpochTDB t, const Eigen::VectorXd& state) {
    if (!std::isfinite(t.mjd())) throw std::runtime_error("Propagator: time is NaN/Inf");
    if (!state.allFinite()) throw std::runtime_error("Propagator: state contains NaN/Inf");

    Eigen::Vector3d position = state.head<3>();
    Eigen::Vector3d velocity = state.tail<3>();
    
    // Use centralized ForceField (on the fly or cached)
    ForceField force(settings_, ephemeris_);
    Eigen::Vector3d acc = force.total_acceleration(t, position, velocity);

    if (!acc.allFinite()) throw std::runtime_error("Propagator: acceleration contains NaN/Inf");

    Eigen::VectorXd xdot(6); xdot << velocity, acc;
    return xdot;
}

Eigen::VectorXd Propagator::integrate_raw_au(const Eigen::VectorXd& y0_au, double t0_mjd, double tf_mjd) {
    DerivativeFunction f = [this](double t_val, const Eigen::VectorXd& y) {
        return compute_derivatives(time::EpochTDB::from_mjd(t_val), y);
    };
    return integrator_->integrate(f, y0_au, t0_mjd, tf_mjd);
}

std::vector<Eigen::VectorXd> Propagator::integrate_raw_au_batch(const Eigen::VectorXd& y0_au, double t0_mjd, const std::vector<double>& tf_mjds) {
    DerivativeFunction f = [this](double t_val, const Eigen::VectorXd& y) {
        return compute_derivatives(time::EpochTDB::from_mjd(t_val), y);
    };
    return integrator_->integrate_at(f, y0_au, t0_mjd, tf_mjds);
}

// Stubs for legacy methods (removed from header?)
void Propagator::update_force_cache(time::EpochTDB t) {}
void Propagator::setup_asteroid_perturbations() {}
Eigen::Vector3d Propagator::compute_n_body_acceleration(const Eigen::Vector3d& position) { return Eigen::Vector3d::Zero(); }
Eigen::Vector3d Propagator::compute_harmonic_acceleration(const Eigen::Vector3d& pos, time::EpochTDB t) { return Eigen::Vector3d::Zero(); }
Eigen::Vector3d Propagator::compute_earth_j2(const Eigen::Vector3d& pos, time::EpochTDB t) { return Eigen::Vector3d::Zero(); }
Eigen::Vector3d Propagator::compute_sun_j2(const Eigen::Vector3d& pos) { return Eigen::Vector3d::Zero(); }
Eigen::Vector3d Propagator::compute_non_gravitational_acceleration(const Eigen::Vector3d& p, const Eigen::Vector3d& v, time::EpochTDB t) { return Eigen::Vector3d::Zero(); }
Eigen::Vector3d Propagator::relativistic_correction(const Eigen::Vector3d& p, const Eigen::Vector3d& v) const { return Eigen::Vector3d::Zero(); }

} // namespace astdyn::propagation
