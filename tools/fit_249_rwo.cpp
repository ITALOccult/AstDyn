/**
 * @file fit_249_rwo.cpp
 * @brief Orbit Determination for Asteroid 249 Ilse using DE441 (Native) and RWO observations
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <cmath>

#include "astdyn/AstDyn.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/ephemeris/PlanetaryData.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/orbit_determination/LeastSquaresFitter.hpp"
#include "astdyn/orbit_determination/ResidualCalculator.hpp"
#include "astdyn/propagation/STMPropagator.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/io/parsers/AstDysRWOParser.hpp"
#include "astdyn/coordinates/KeplerianElements.hpp"
#include "astdyn/coordinates/CartesianState.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/io/parsers/OrbFitEQ1Parser.hpp"

using namespace astdyn;
using namespace astdyn::orbit_determination;
using namespace astdyn::propagation;
using namespace astdyn::ephemeris;
using namespace astdyn::coordinates; 

// Configuration
const std::string RWO_FILE = "249.rwo";
const std::string SPK_FILE = "/Users/michelebigi/Downloads/de441_part-2.bsp";
const bool USE_RELATIVITY = true;
const bool USE_ASTEROIDS = true; 

// Units Constants
const double AU_IN_KM = 149597870.7;

int main() {
    try {
        std::cout << "=== Asteroid 249 Ilse Orbit Determination ===\n\n";

        // 1. Initialize Ephemeris Provider (DE441)
        std::cout << "Loading DE441 Ephemeris (" << SPK_FILE << ")... ";
        auto de441 = std::make_shared<DE441Provider>(SPK_FILE);
        PlanetaryEphemeris::setProvider(de441); 
        auto ephem_ptr = std::make_shared<PlanetaryEphemeris>(); 
        std::cout << "DONE.\n";

        // 2. Load Observations
        std::cout << "Loading RWO Observations (" << RWO_FILE << ")... ";
        io::parsers::AstDysRWOParser parser;
        auto observations = parser.parse(RWO_FILE);
        std::cout << "DONE. Loaded " << observations.size() << " observations.\n";
        
        using OpticalObs = astdyn::io::IObservationParser::OpticalObservation;
        
        // Filter observations near initial epoch (2017) to ensure convergence
        std::vector<OpticalObs> valid_obs;
        // MJD 58016 is Sept 2017. Let's take +/- 100 days.
        double min_mjd = 57900.0; 
        double max_mjd = 58150.0;
        
        for (const auto& o : observations) {
            if (o.mjd_utc >= min_mjd && o.mjd_utc <= max_mjd) {
                valid_obs.push_back(o);
            }
        }
        
        if (valid_obs.empty()) {
            std::cerr << "Warning: No obs in 2017 range. Using all post-2010.\n";
            // Fallback
             min_mjd = 55197.0; 
             for (const auto& o : observations) { if (o.mjd_utc >= min_mjd) valid_obs.push_back(o); }
        }

        // Use subset of observations (all found in range, up to 100)
        std::vector<OpticalObs> obs_subset;
        size_t n_obs = 100;
        if (valid_obs.size() > n_obs) {
             // Take evenly spaced or first N? Let's take first N since close to epoch
             obs_subset.insert(obs_subset.end(), valid_obs.begin(), valid_obs.begin() + n_obs);
        } else {
             obs_subset = valid_obs;
        }
        
        std::cout << "Using subset of " << obs_subset.size() << " observations from 2017 window.\n";
        
        // Disable epoch shift since we fit AT epoch 58016
        bool shift_epoch = false; 
        
        // Variable "fit_epoch" will be handling this below?
        // We must ensure fit logic uses t0 if no shift.
        // Looking at code below:
        // if (shift_epoch) ... 
        // else x0 = cart_au...
        // This is handled logic.

        // 3. Initial Guess from 249.eq1 (AstDyS Ecliptic Elements)
        std::cout << "Loading Initial Elements from 249.eq1...\n";
        astdyn::io::parsers::OrbFitEQ1Parser eq1_parser;
        // Assume 249.eq1 is in current directory tools/
        auto initial_elem = eq1_parser.parse("249.eq1");
        
        // Convert OrbitalElements (struct) to KeplerianElements (class)
        // Adjust if types differ. OrbitalElements struct has fields similar to KeplerianElements.
        // We need astdyn::coordinates::KeplerianElements class for conversion.
        astdyn::coordinates::KeplerianElements initial_kep;
        initial_kep.set_semi_major_axis(initial_elem.semi_major_axis * AU_IN_KM); 
        initial_kep.set_eccentricity(initial_elem.eccentricity);
        initial_kep.set_inclination(initial_elem.inclination); 
        initial_kep.set_RAAN(initial_elem.longitude_asc_node); 
        initial_kep.set_argument_of_periapsis(initial_elem.argument_perihelion); 
        initial_kep.set_mean_anomaly(initial_elem.mean_anomaly); 
        
        double t0 = initial_elem.epoch_mjd_tdb;
        std::cout << "Loaded Elements @ MJD " << t0 << " (Ecliptic J2000)\n";
        
        // 4. Setup Propagator Settings (moved up to set GM for KeplerianElements)
        astdyn::propagation::PropagatorSettings settings; 
        // Assuming Sun's GM (km^3/s^2) for initial Keplerian elements
        settings.central_body_gm = 1.327124400419394e11; 
        initial_kep.set_mu(settings.central_body_gm);

        // Convert to Cartesian (Ecliptic)
        CartesianState cart_km_ecl = initial_kep.to_cartesian(); 
        
        // ROTATION: Ecliptic -> Equatorial J2000 (ICRF)
        // AstDyS .eq1 header says "ECLM J2000", Propagator needs Equatorial.
        double eps = 23.4392911 * M_PI / 180.0; // Mean Obliquity J2000
        double cos_eps = std::cos(eps);
        double sin_eps = std::sin(eps);
        
        auto rotate_ecl_to_eq = [&](const Eigen::Vector3d& v_ecl) {
            return Eigen::Vector3d(
                v_ecl.x(),
                v_ecl.y() * cos_eps - v_ecl.z() * sin_eps,
                v_ecl.y() * sin_eps + v_ecl.z() * cos_eps
            );
        };
        
        CartesianState cart_km_eq;
        cart_km_eq.set_position(rotate_ecl_to_eq(cart_km_ecl.position()));
        cart_km_eq.set_velocity(rotate_ecl_to_eq(cart_km_ecl.velocity()));
        
        // Convert to AU/day
        CartesianState cart_au = cart_km_eq.to_AU_per_day();
        
        // Initial for fitter
        Eigen::Vector<double, 6> x0;
        x0.head<3>() = cart_au.position();
        x0.tail<3>() = cart_au.velocity();

        std::cout << "Initial State (AU, AU/d) @ MJD " << t0 << " (Rotated to EQ): \n" << x0.transpose() << "\n\n";

        // 4. Setup Propagator Settings (rest of settings)
        // CRITICAL: Propagator works in AU and Days, so GM must be in AU^3/day^2.
        // We used km/s for Keplerian conversion, but now switch to AU/day.
        settings.central_body_gm = 0.0002959122082855911; // GM Sun in AU^3/d^2
        
        // Full High-Precision Model
        settings.include_planets = true;
        settings.perturb_mercury = true;
        settings.perturb_venus = true;
        settings.perturb_earth = true;
        settings.perturb_mars = true;
        settings.perturb_jupiter = true;
        settings.perturb_saturn = true;
        settings.perturb_uranus = true;
        settings.perturb_neptune = true;
        settings.include_relativity = USE_RELATIVITY;
        settings.include_asteroids = USE_ASTEROIDS;

        // Create Propagator for setup/shift
        auto createPropagator = [&]() {
            auto integ = std::make_unique<RKF78Integrator>(0.1, 1e-12);
            return std::make_shared<Propagator>(std::move(integ), ephem_ptr, settings);
        };
        auto setup_prop = createPropagator();

        // --- EPOCH SHIFT STRATEGY ---
        // Propagate state from 2017 (t0) to the first observation time (2024).
        // This reduces linearization error and prevents divergence.
        double t_first_obs = obs_subset.front().mjd_utc;
        std::cout << "\n[Step 3b] Propagating initial guess to new epoch: " << t_first_obs << " (MJD)...\n";
        
        // Convert Vector6d (x0) -> CartesianElements -> Propagate -> Vector6d (x0_new)
        // Note: x0 is in AU, AU/d. t0 is 58016.
        CartesianState cs_t0(x0.head<3>(), x0.tail<3>(), cart_au.mu());
        
        // Prepare propagation
        CartesianElements el_in;
        el_in.position = cs_t0.position();
        el_in.velocity = cs_t0.velocity();
        el_in.epoch_mjd_tdb = t0;
        el_in.gravitational_parameter = settings.central_body_gm; // FIX: Prevent Zero Acc

        // Propagate
        CartesianElements el_out = setup_prop->propagate_cartesian(el_in, t_first_obs);
        
        // Update x0 and t0 for the fitter
        x0.head<3>() = el_out.position;
        x0.tail<3>() = el_out.velocity;
        t0 = t_first_obs;

        std::cout << "New Initial State (AU, AU/d) @ Epoch " << t0 << ":\n" << x0.transpose() << "\n\n";

        // 4. Setup STM Force Model (reusing master_prop created below or creating new one)
        auto master_prop = createPropagator();
        
        // DEBUG: Check Physics
        std::cout << "DEBUG: Central Body GM = " << settings.central_body_gm << "\n";
        auto acc_check = master_prop->compute_derivatives(t0, x0);
        std::cout << "DEBUG: Initial Acc @ t0: " << acc_check.tail<3>().transpose() << " | Norm: " << acc_check.tail<3>().norm() << "\n";

        STMPropagator::ForceFunction force_func = [master_prop](double t, const Eigen::Vector<double, 6>& state) {
            return master_prop->compute_derivatives(t, state); 
        };

        // 5. Initialize Fitter
        LeastSquaresFitter fitter;
        fitter.set_tolerance(1e-6);
        fitter.set_max_iterations(5);
        fitter.set_outlier_threshold(3.0); 

        // Residual Function
        LeastSquaresFitter::ResidualFunction res_func = [&](const Eigen::Vector<double, 6>& state_at_epoch, double epoch) {
            std::vector<ObservationResidual> residuals_vec;
            residuals_vec.reserve(obs_subset.size());
            
            auto prop = createPropagator();
            
            CartesianState curr_cs(state_at_epoch.head<3>(), state_at_epoch.tail<3>(), cart_au.mu());
            double curr_time = epoch;
            
            ResidualCalculator res_calculator;
            res_calculator.set_ephemeris_provider(de441);
            res_calculator.load_observatories("./observatories.dat");
            std::cout << "DEBUG: Loaded observatories from ./observatories.dat\n";
            
            for (const auto& obs : obs_subset) {
                double t_target = obs.mjd_utc + (37.0 + 32.184)/86400.0;
                
                // Debug Topocentric Link
                static bool first_debug = true;
                if (first_debug) {
                    // std::cout << "DEBUG FRAME: Processing Obs Code " << obs.obs_code << "\n";
                    // Check if calculator has it (implied by computation result later)
                    first_debug = false;
                }
                
                if (std::abs(t_target - curr_time) > 1e-9) {
                    CartesianElements el_in;
                    el_in.position = curr_cs.position();
                    el_in.velocity = curr_cs.velocity();
                    el_in.epoch_mjd_tdb = curr_time; // Fix: Set epoch to prevent 1858 default
                    el_in.gravitational_parameter = cart_au.mu(); // FIX: Prevent Zero Acc
                    
                    auto el_out = prop->propagate_cartesian(el_in, t_target);
                    
                    curr_cs.set_position(el_out.position);
                    curr_cs.set_velocity(el_out.velocity);
                    curr_time = t_target;
                }
                
                CartesianState curr_km = curr_cs.to_km_per_s();
                Eigen::Vector<double, 6> obj_km;
                obj_km.head<3>() = curr_km.position();
                obj_km.tail<3>() = curr_km.velocity();
                
                astdyn::orbit_determination::Observation fit_obs;
                fit_obs.ra_deg = obs.ra * 180.0/M_PI; 
                fit_obs.dec_deg = obs.dec * 180.0/M_PI;
                fit_obs.epoch_mjd = obs.mjd_utc;
                fit_obs.observatory_code = obs.obs_code;
                fit_obs.ra_sigma_arcsec = (obs.sigma_ra > 0) ? obs.sigma_ra : 1.0;
                fit_obs.dec_sigma_arcsec = (obs.sigma_dec > 0) ? obs.sigma_dec : 1.0;
                fit_obs.weight = 1.0; 
                fit_obs.rejected = false;

                auto res_base = res_calculator.compute_residual(fit_obs, obj_km, t_target); 
                
                ObservationResidual full_res;
                full_res.epoch_mjd = res_base.epoch_mjd;
                full_res.ra_computed_deg = res_base.ra_computed_deg;
                full_res.dec_computed_deg = res_base.dec_computed_deg;
                full_res.ra_residual_arcsec = res_base.ra_residual_arcsec;
                full_res.dec_residual_arcsec = res_base.dec_residual_arcsec;
                full_res.rejected = res_base.rejected;
                full_res.ra_obs_deg = fit_obs.ra_deg;
                full_res.dec_obs_deg = fit_obs.dec_deg;
                full_res.weight_ra = 1.0;
                full_res.weight_dec = 1.0;

                residuals_vec.push_back(full_res);
            }
            return residuals_vec;
        };

        // STM Function
        LeastSquaresFitter::STMFunction stm_func = [&](const Eigen::Vector<double, 6>& state_at_epoch, double t0_func, double tf_func) {
            auto integ = std::make_unique<RKF78Integrator>(0.1, 1e-9);
            STMPropagator stm_prop(std::move(integ), force_func, nullptr);
            
            Eigen::Matrix<double, 6, 6> Phi0 = Eigen::Matrix<double, 6, 6>::Identity();
            auto result = stm_prop.propagate(state_at_epoch, t0_func, tf_func, Phi0);
            return std::make_pair(result.state, result.stm);
        };

        std::cout << "Starting Least Squares Fit...\n";
        auto result = fitter.fit(x0, t0, res_func, stm_func);

        std::cout << "\n=== Fit Results ===\n";
        std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
        std::cout << "RMS RA:  " << result.rms_ra_arcsec << " arcsec\n";
        std::cout << "RMS Dec: " << result.rms_dec_arcsec << " arcsec\n";
        std::cout << "RMS Tot: " << result.rms_total_arcsec << " arcsec\n";
        std::cout << "\nFitted State (AU, AU/d) @ MJD " << t0 << ":\n" << result.state.transpose() << "\n";
        
        CartesianState fit_au(result.state.head<3>(), result.state.tail<3>(), cart_au.mu());
        CartesianState fit_km = fit_au.to_km_per_s();
        auto final_kep = astdyn::coordinates::KeplerianElements::from_cartesian(fit_km);

        std::cout << "\nFitted Elements:\n";
        std::cout << "  a: " << final_kep.semi_major_axis() / AU_IN_KM << " AU\n";
        std::cout << "  e: " << final_kep.eccentricity() << "\n";
        std::cout << "  i: " << final_kep.inclination() * 180.0/M_PI << " deg\n";
        std::cout << "  Node: " << final_kep.RAAN() * 180.0/M_PI << " deg\n";
        std::cout << "  Peri: " << final_kep.argument_of_periapsis() * 180.0/M_PI << " deg\n";
        std::cout << "  Mean Anom: " << final_kep.mean_anomaly() * 180.0/M_PI << " deg\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
