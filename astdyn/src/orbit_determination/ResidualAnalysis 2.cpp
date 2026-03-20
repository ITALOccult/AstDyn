/**
 * @file ResidualAnalysis.cpp
 * @brief Analysis of orbit residuals (O-C).
 */

#include "astdyn/orbit_determination/ResidualAnalysis.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/time/TimeScale.hpp"
#include <iomanip>
#include <sstream>
#include <cmath>

namespace astdyn::orbit_determination {

OrbitValidationSummary ResidualAnalysis::analyze_orbit(
    const physics::CartesianStateTyped<core::GCRF>& state,
    const std::vector<observations::OpticalObservation>& obs,
    std::shared_ptr<propagation::Propagator> propagator)
{
    OrbitValidationSummary summary;
    std::stringstream ss;

    ss << "--- Residual Analysis Report ---\n";
    ss << std::fixed << std::setprecision(3);
    ss << "Observation          Time (TDB)    dRA (\")    dDec (\")\n";
    ss << "---------------------------------------------------------\n";

    if (obs.empty()) {
        summary.rms_ra = 0.0;
        summary.rms_dec = 0.0;
        summary.rms_total = 0.0;
        summary.report_text = "No observations to analyze.";
        return summary;
    }

    double sum_sq_ra = 0.0;
    double sum_sq_dec = 0.0;

    auto residual_calc = std::make_unique<ResidualCalculator<core::GCRF>>(propagator->get_ephemeris(), propagator);

    for (size_t i = 0; i < obs.size(); ++i) {
        const auto& o = obs[i];
        time::EpochTDB t_obs = time::to_tdb(o.time);

        // Propagate state to observation time
        auto state_at_t = propagator->propagate_cartesian(state, t_obs);

        // Get proper observer position (including station offset)
        auto obs_pos_opt = residual_calc->get_observer_position(o);
        if (!obs_pos_opt) continue;
        
        // Calculate Line of Sight (relative to observer)
        Eigen::Vector3d rho_vec = state_at_t.position.to_eigen_si() - obs_pos_opt->to_eigen_si();

        double rho = rho_vec.norm();
        double ra_calc = std::atan2(rho_vec.y(), rho_vec.x());
        double dec_calc = std::asin(std::max(-1.0, std::min(1.0, rho_vec.z() / rho)));

        // Residuals
        double d_ra = o.ra.to_rad() - ra_calc;
        while (d_ra > constants::PI) d_ra -= constants::TWO_PI;
        while (d_ra < -constants::PI) d_ra += constants::TWO_PI;
        double d_dec = o.dec.to_rad() - dec_calc;

        // Convert to arcseconds (Projected RA)
        double dra_arcsec = d_ra * std::cos(dec_calc) * constants::RAD_TO_ARCSEC;
        double ddec_arcsec = d_dec * constants::RAD_TO_ARCSEC;

        ss << std::setw(3) << i << "  " << t_obs.mjd() << "  " 
           << std::setw(8) << dra_arcsec << "  " << std::setw(8) << ddec_arcsec << "\n";

        sum_sq_ra += dra_arcsec * dra_arcsec;
        sum_sq_dec += ddec_arcsec * ddec_arcsec;
    }

    summary.rms_ra = std::sqrt(sum_sq_ra / obs.size());
    summary.rms_dec = std::sqrt(sum_sq_dec / obs.size());
    summary.rms_total = std::sqrt((sum_sq_ra + sum_sq_dec) / (2.0 * obs.size()));

    ss << "---------------------------------------------------------\n";
    ss << "RMS RA (proj):  " << summary.rms_ra << " \"\n";
    ss << "RMS Dec:        " << summary.rms_dec << " \"\n";
    ss << "RMS Total:      " << summary.rms_total << " \"\n";

    summary.report_text = ss.str();
    return summary;
}

} // namespace astdyn::orbit_determination
