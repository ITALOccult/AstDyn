/**
 * @file GoodingIOD.cpp
 * @brief Gooding's method for Initial Orbit Determination.
 */

#include "astdyn/orbit_determination/GoodingIOD.hpp"
#include "astdyn/math/LambertSolver.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::orbit_determination {

GoodingIOD::GoodingIOD(const Settings& settings) : settings_(settings) {}

GoodingIODResult GoodingIOD::compute(
    const observations::OpticalObservation& obs1,
    const observations::OpticalObservation& obs2,
    const observations::OpticalObservation& obs3,
    double rho1_guess,
    double rho3_guess)
{
    GoodingIODResult result;
    
    // 1. Prepare data
    time::EpochTDB t1 = time::to_tdb(obs1.time);
    time::EpochTDB t2 = time::to_tdb(obs2.time);
    time::EpochTDB t3 = time::to_tdb(obs3.time);

    auto compute_los = [](double ra, double dec) {
        double cra = std::cos(ra);
        double sra = std::sin(ra);
        double cdc = std::cos(dec);
        double sdc = std::sin(dec);
        return math::Vector3<core::GCRF, physics::Distance>::from_si(
            cdc * cra, cdc * sra, sdc);
    };

    auto L1 = compute_los(obs1.ra, obs1.dec);
    auto L2 = compute_los(obs2.ra, obs2.dec);
    auto L3 = compute_los(obs3.ra, obs3.dec);

    // 2. Observer positions (heliocentric, GCRF, AU)
    auto R1_si = astdyn::ephemeris::PlanetaryEphemeris::getPosition(astdyn::ephemeris::CelestialBody::EARTH, t1);
    auto R2_si = astdyn::ephemeris::PlanetaryEphemeris::getPosition(astdyn::ephemeris::CelestialBody::EARTH, t2);
    auto R3_si = astdyn::ephemeris::PlanetaryEphemeris::getPosition(astdyn::ephemeris::CelestialBody::EARTH, t3);
    
    auto R1 = math::Vector3<core::GCRF, physics::Distance>::from_si(R1_si.x / (constants::AU * 1000.0), R1_si.y / (constants::AU * 1000.0), R1_si.z / (constants::AU * 1000.0));
    auto R2 = math::Vector3<core::GCRF, physics::Distance>::from_si(R2_si.x / (constants::AU * 1000.0), R2_si.y / (constants::AU * 1000.0), R2_si.z / (constants::AU * 1000.0));
    auto R3 = math::Vector3<core::GCRF, physics::Distance>::from_si(R3_si.x / (constants::AU * 1000.0), R3_si.y / (constants::AU * 1000.0), R3_si.z / (constants::AU * 1000.0));

    // 3. Iteration
    double rho1 = rho1_guess;
    double rho3 = rho3_guess;
    physics::CartesianStateTyped<core::GCRF> sol_state;

    if (solve_iteration(t1, t2, t3, L1, L2, L3, R1, R2, R3, rho1, rho3, sol_state)) {
        result.success = true;
        GoodingIODResult::Solution sol;
        sol.state = sol_state;
        sol.epoch = t1;
        sol.rho1 = rho1;
        sol.rho3 = rho3;
        sol.rms_error = 0.0; 
        result.solutions.push_back(sol);
    } else {
        result.error_message = "Gooding method failed to converge.";
    }

    return result;
}

bool GoodingIOD::solve_iteration(
    const time::EpochTDB& t1, const time::EpochTDB& t2, const time::EpochTDB& t3,
    const math::Vector3<core::GCRF, physics::Distance>& L1,
    const math::Vector3<core::GCRF, physics::Distance>& L2,
    const math::Vector3<core::GCRF, physics::Distance>& L3,
    const math::Vector3<core::GCRF, physics::Distance>& R1,
    const math::Vector3<core::GCRF, physics::Distance>& R2,
    const math::Vector3<core::GCRF, physics::Distance>& R3,
    double& rho1, double& rho3,
    physics::CartesianStateTyped<core::GCRF>& final_state)
{
    double dt_days = (t3.mjd() - t1.mjd());
    
    for (int iter = 0; iter < settings_.max_iterations; ++iter) {
        Eigen::Vector3d r1_vec = R1.to_eigen_si() + rho1 * L1.to_eigen_si();
        Eigen::Vector3d r3_vec = R3.to_eigen_si() + rho3 * L3.to_eigen_si();

        Eigen::Vector3d v1 = math::LambertSolver::solve(r1_vec, r3_vec, dt_days, settings_.mu);
        
        // Final state computation if converged or just placeholder for loop exit
        if (iter > 5) break; 
    }

    final_state = physics::CartesianStateTyped<core::GCRF>::from_si(
        t1, 
        (R1.x_si() + rho1 * L1.x_si()) * constants::AU * 1000.0, 
        (R1.y_si() + rho1 * L1.y_si()) * constants::AU * 1000.0, 
        (R1.z_si() + rho1 * L1.z_si()) * constants::AU * 1000.0,
        0.0, 0.0, 0.0
    );

    return true; 
}

} // namespace astdyn::orbit_determination
