/**
 * @file ResidualAnalysis.cpp
 * @brief Implementation of orbit validation summary report.
 */

#include "astdyn/orbit_determination/ResidualAnalysis.hpp"
#include <iomanip>
#include <sstream>
#include <cmath>

namespace astdyn::orbit_determination {

OrbitValidationSummary ResidualAnalysis::analyze_orbit(
    const physics::CartesianStateTyped<core::GCRF>& initial_state,
    const std::vector<observations::OpticalObservation>& obs,
    std::shared_ptr<propagation::Propagator> propagator)
{
    OrbitValidationSummary results;
    results.num_observations = (int)obs.size();
    
    std::stringstream report;
    report << "--- Orbit Validation Report ---\n";
    report << std::left << std::setw(20) << "Time (MJD)" 
           << std::setw(15) << "RA-res (\")" 
           << std::setw(15) << "Dec-res (\")" << "\n";
    report << std::string(50, '-') << "\n";

    double sum_sq_ra = 0.0, sum_sq_dec = 0.0, sum_sq_tot = 0.0;
    
    for (const auto& o : obs) {
        time::EpochTDB t_obs = time::to_tdb(o.time);
        
        // Propagate
        auto state_calc = propagator->propagate_cartesian(initial_state, t_obs);
        
        // Predict observation (geocentric)
        auto earth_pos_si = ephemeris::PlanetaryEphemeris::getPosition(ephemeris::CelestialBody::EARTH, t_obs);
        auto R_obs = earth_pos_si.to_eigen(); 
        
        // rho = r - R
        auto rho_vec = state_calc.position.to_eigen_si() - R_obs;
        double rho = rho_vec.norm();
        double ra_calc = std::atan2(rho_vec.y(), rho_vec.x());
        double dec_calc = std::asin(rho_vec.z() / rho);
        
        // Differences (arcsec)
        double dra = (o.ra - ra_calc);
        while (dra > M_PI) dra -= 2.0 * M_PI;
        while (dra < -M_PI) dra += 2.0 * M_PI;
        dra *= std::cos(o.dec) * constants::RAD_TO_DEG * 3600.0;
        
        double ddec = (o.dec - dec_calc) * constants::RAD_TO_DEG * 3600.0;
        
        sum_sq_ra += dra * dra;
        sum_sq_dec += ddec * ddec;
        sum_sq_tot += dra * dra + ddec * ddec;
        
        report << std::fixed << std::setprecision(5) << std::setw(20) << o.time.mjd() 
               << std::setw(15) << dra 
               << std::setw(15) << ddec << "\n";
    }

    results.rms_ra = std::sqrt(sum_sq_ra / results.num_observations);
    results.rms_dec = std::sqrt(sum_sq_dec / results.num_observations);
    results.rms_total = std::sqrt(sum_sq_tot / results.num_observations);

    report << std::string(50, '-') << "\n";
    report << "RMS RA:   " << std::fixed << std::setprecision(3) << results.rms_ra << " arcsec\n";
    report << "RMS Dec:  " << std::fixed << std::setprecision(3) << results.rms_dec << " arcsec\n";
    report << "RMS Total:" << std::fixed << std::setprecision(3) << results.rms_total << " arcsec\n";
    
    results.report_text = report.str();
    return results;
}

} // namespace astdyn::orbit_determination
