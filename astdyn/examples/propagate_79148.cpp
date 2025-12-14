/**
 * @file propagate_79148.cpp
 * @brief Propagate Asteroid 79148 to 2026-04-30 for JPL verification.
 * 
 * Features:
 * - Load initial state from 79148.eq1
 * - Full Physics: Sun, 8 Planets, Relativity (PPN)
 * - Target Epoch: 2026-04-30 TDB
 * - Output: Cartesian State (Equatorial J2000) and Keplerian Elements.
 */

#include <astdyn/AstDynEngine.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/time/TimeScale.hpp>
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <astdyn/coordinates/CartesianState.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace astdyn;

// Helper to trim leading whitespace (same as in simple_pompeja_fit)
std::string ltrim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t");
    return (start == std::string::npos) ? "" : s.substr(start);
}

// Parse OrbFit .eq1 file into propagation::EquinoctialElements
propagation::EquinoctialElements parse_eq1_file(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    propagation::EquinoctialElements equ;
    std::string line;
    bool found_equ = false;
    bool found_mjd = false;

    while (std::getline(file, line)) {
        std::string trimmed = ltrim(line);
        if (trimmed.empty() || trimmed[0] == '!') continue;

        if (trimmed.substr(0, 3) == "EQU") {
            std::istringstream iss(trimmed.substr(3));
            iss >> equ.a >> equ.h >> equ.k >> equ.p >> equ.q >> equ.lambda;
            // Convert mean longitude from degrees to radians
            equ.lambda *= constants::PI / 180.0;
            found_equ = true;
        }
        else if (trimmed.substr(0, 3) == "MJD") {
            std::istringstream iss(trimmed.substr(3));
            iss >> equ.epoch_mjd_tdb;
            found_mjd = true;
        }

        if (found_equ && found_mjd) break;
    }

    if (!found_equ || !found_mjd) {
        throw std::runtime_error("Could not parse equinoctial elements from file");
    }

    std::cerr << "Parsed Raw Elements:\n";
    std::cerr << "  a = " << equ.a << "\n";
    std::cerr << "  h = " << equ.h << "\n";
    std::cerr << "  k = " << equ.k << "\n";
    std::cerr << "  p = " << equ.p << "\n";
    std::cerr << "  q = " << equ.q << "\n";
    std::cerr << "  lambda = " << equ.lambda * 180.0 / constants::PI << " deg\n";

    equ.gravitational_parameter = constants::GMS;
    return equ;
}

int main(int argc, char** argv) {
    try {
        std::cout << "=== Asteroid 79148 Propagation Verification (User Elements) ===\n";
        
        // 1. User Provided Elements (Keplerian Ecliptic J2000)
        // EPOCH=  2458749.5 ! 2019-Sep-23.00 (TDB)
        // EC= .1693869472540817   QR= 2.31958229445229    TP= 2459553.8582709809      
        // OM= 148.2842362019823   W=  312.8020637281074   IN= 4.750216031450175 
        
        // 1. Initial State provided by OrbFit .eq1 or .ele file (example placeholder)
        // Here we use the values for 79148 again but cleaner.
        // In a real example, this might parse argv[1].
        
        // Example: Asteroid 79148 (1992 SN3)
        // Epoch: 58749.5 MJD TDB
        propagation::CartesianElements cart_eq;
        cart_eq.epoch_mjd_tdb = time::jd_to_mjd(2458749.5);
        cart_eq.gravitational_parameter = constants::GMS;
        
        // Horizons State Vector (Heliocentric ICRF)
        cart_eq.position = Eigen::Vector3d(1.032300572794974, -2.901750550202138, -1.069323490646001);
        cart_eq.velocity = Eigen::Vector3d(0.008162881713437332, 0.002917128256631872, 0.0006499392031295894);
        
        // For display
        auto kep_eq = propagation::cartesian_to_keplerian(cart_eq);
        
        std::cout << "Initial Epoch: " << std::fixed << std::setprecision(5) << kep_eq.epoch_mjd_tdb << " MJD TDB\n";
        std::cout << "Initial State (Equatorial J2000):\n";
        std::cout << "  X = " << cart_eq.position[0] << " AU\n";
        std::cout << "  Y = " << cart_eq.position[1] << " AU\n";
        std::cout << "  Z = " << cart_eq.position[2] << " AU\n";
        std::cout << " VX = " << cart_eq.velocity[0] << " AU/d\n";
        std::cout << " VY = " << cart_eq.velocity[1] << " AU/d\n";
        std::cout << " VZ = " << cart_eq.velocity[2] << " AU/d\n\n";

        // 2. Setup Engine
        AstDynEngine engine;
        AstDynConfig config = engine.config();
        
        // HIGH FIDELITY SETTINGS
        config.integrator_type = "RKF78";
        config.initial_step_size = 0.5;
        config.tolerance = 1e-13;
        
        config.propagator_settings.include_planets = true;
        config.propagator_settings.perturb_mercury = true;
        config.propagator_settings.perturb_venus = true;
        config.propagator_settings.perturb_earth = true;
        config.propagator_settings.perturb_mars = true;
        config.propagator_settings.perturb_jupiter = true;
        config.propagator_settings.perturb_saturn = true;
        config.propagator_settings.include_relativity = true;
        config.propagator_settings.include_asteroids = true;
        
        config.ephemeris_type = "DE441"; 
        config.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
        config.asteroid_ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/codes_300ast_20100725.bsp";
        
        engine.set_config(config);
        engine.set_initial_orbit(kep_eq);

        // 3. Propagate to Target Date
        // Target: 2026-04-30 00:00:00 TDB
        double target_mjd = time::calendar_to_mjd(2026, 4, 30, 0.0);
        
        std::cout << "Propagating to " << target_mjd << " MJD (2026-04-30 00:00 UTC/TDB)...\n";
        
        auto final_state = engine.propagate_to(target_mjd);
        auto final_cart = propagation::keplerian_to_cartesian(final_state);
        
        std::cout << "\n=== FINAL STATE (2026-04-30 00:00:00) ===\n";
        std::cout << "Epoch: " << final_cart.epoch_mjd_tdb << " MJD TDB\n";
        
        std::cout << std::setprecision(12);
        std::cout << "Position (AU) [J2000 Eq]:\n";
        std::cout << " X = " << final_cart.position[0] << "\n";
        std::cout << " Y = " << final_cart.position[1] << "\n";
        std::cout << " Z = " << final_cart.position[2] << "\n";
        
        std::cout << "Velocity (AU/d) [J2000 Eq]:\n";
        std::cout << " VX = " << final_cart.velocity[0] << "\n";
        std::cout << " VY = " << final_cart.velocity[1] << "\n";
        std::cout << " VZ = " << final_cart.velocity[2] << "\n";
        
        std::cout << "\nKeplerian Elements (J2000 Eq):\n";
        std::cout << " a = " << final_state.semi_major_axis << " AU\n";
        std::cout << " e = " << final_state.eccentricity << "\n";
        std::cout << " i = " << final_state.inclination * constants::RAD_TO_DEG << " deg\n";
        std::cout << " Ω = " << final_state.longitude_ascending_node * constants::RAD_TO_DEG << " deg\n";
        std::cout << " ω = " << final_state.argument_perihelion * constants::RAD_TO_DEG << " deg\n";
        std::cout << " M = " << final_state.mean_anomaly * constants::RAD_TO_DEG << " deg\n";

        // 4. Compute Geocentric Equatorial Coordinates
        std::cout << "\n=== GEOCENTRIC COORDINATES (Astrometric J2000) ===\n";
        
        // Get Earth Position at Target Epoch
        // Note: AstDynEngine initializes Ephemeris. 
        // DE441Provider (Native) returns Equatorial J2000 directly.
        // PlanetaryEphemeris::getState forwards this.
        // Analytical returns Ecliptic.
        // Since we are forcing DE441 now, we assume Equatorial output.
        
        double target_jd_tdb = time::mjd_to_jd(target_mjd);
        auto earth_state_raw = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, target_jd_tdb);
        
        // Check if we need to transform (Analytical = Ecliptic, DE441 = Equatorial)
        // For this specific verification with DE441 loaded:
        Vector3d earth_pos_bary = earth_state_raw.position();
        
        // CORRECTION: Convert Earth Barycentric -> Heliocentric
        // Propagator integrates in Heliocentric frame. DE441 returns Barycentric.
        // Earth_Helio = Earth_Bary - Sun_Bary
        auto sun_state_raw = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::SUN, target_jd_tdb);
        Vector3d sun_pos_bary = sun_state_raw.position();
        
        Vector3d earth_pos = earth_pos_bary - sun_pos_bary;
        
        // If we were using Analytical, we would need:
        // auto earth_eq = coordinates::ReferenceFrame::transform_state(earth_state_raw, ECLIPTIC, J2000);
        // earth_pos = earth_eq.position();
        
        // Correzione per il Tempo Luce (Light-Time Correction)
        // JPL "Astrometric" include solitamente la correzione per il tempo luce.
        // Dobbiamo trovare la posizione dell'asteroide al tempo (t - tau)
        
        double c_au_day = 173.144632674240; // Velocità della luce in AU/d
        double tau = 0.0;
        Vector3d asteroid_pos_retarded = final_cart.position;
        double range = (asteroid_pos_retarded - earth_pos).norm();
        
        // Iterazione per il tempo luce (3 iterazioni sono sufficienti)
        for(int k=0; k<3; ++k) {
            tau = range / c_au_day;
            double retarded_mjd = target_mjd - tau;
            
            // Propaghiamo l'asteroide al tempo ritardato
            auto state_retarded = engine.propagate_to(retarded_mjd);
            auto cart_retarded = propagation::keplerian_to_cartesian(state_retarded);
            asteroid_pos_retarded = cart_retarded.position;
            
            range = (asteroid_pos_retarded - earth_pos).norm();
        }
        
        Vector3d rho_vec = asteroid_pos_retarded - earth_pos;
        
        double ra_rad = std::atan2(rho_vec[1], rho_vec[0]);
        if (ra_rad < 0) ra_rad += 2.0 * constants::PI;
        double dec_rad = std::asin(rho_vec[2] / range);
        
        double ra_deg = ra_rad * constants::RAD_TO_DEG;
        double dec_deg = dec_rad * constants::RAD_TO_DEG;
        
        // Format to HMS / DMS
        int ra_h = static_cast<int>(ra_deg / 15.0);
        int ra_m = static_cast<int>((ra_deg / 15.0 - ra_h) * 60.0);
        double ra_s = ((ra_deg / 15.0 - ra_h) * 60.0 - ra_m) * 60.0;
        
        int dec_d = static_cast<int>(std::abs(dec_deg));
        int dec_m = static_cast<int>((std::abs(dec_deg) - dec_d) * 60.0);
        double dec_s = ((std::abs(dec_deg) - dec_d) * 60.0 - dec_m) * 60.0;
        if (dec_deg < 0) dec_d = -dec_d;

        std::cout << "Terra (Elio):      " << earth_pos.transpose() << " AU\n";
        std::cout << "Asteroide (t-tau): " << asteroid_pos_retarded.transpose() << " AU\n";
        std::cout << "Tempo Luce:        " << (tau * 1440.0) << " minuti\n";
        std::cout << "Range:             " << range << " AU\n";
        std::cout << "RA:                " << ra_deg << " deg\n";
        std::cout << "Dec:               " << dec_deg << " deg\n";
        std::cout << "RA (HMS):          " << ra_h << "h " << ra_m << "m " << std::fixed << std::setprecision(3) << ra_s << "s\n";
        std::cout << "Dec (DMS):         " << (dec_deg < 0 && dec_d == 0 ? "-" : "") << dec_d << "d " << dec_m << "m " << dec_s << "s\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
