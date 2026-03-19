#include "astdyn/AstDyn.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace std;

int main() {
    try {
        io::HorizonsClient horizons;
        auto t_obs = time::EpochTDB::from_jd(2461164.5); // May 4, 2026
        
        // 1. Query JPL OBSERVER (Quantity 1: Astrometric RA/Dec)
        cout << "Querying JPL OBSERVER (Quantity 1)..." << endl;
        auto obs_opt = horizons.query_observation("234", t_obs, "500");
        if (!obs_opt) return 1;
        auto obs_jpl = obs_opt.value();
        
        // 2. Query JPL BARYCENTRIC VECTORS (Quantity 2)
        cout << "Querying JPL BARYCENTRIC VECTORS (Quantity 2)..." << endl;
        auto ast_opt = horizons.query_vectors("234", t_obs, "@0");
        auto earth_opt = horizons.query_vectors("399", t_obs, "@0");
        if (!ast_opt || !earth_opt) return 1;
        
        auto p_ast = ast_opt.value().position.to_eigen_si();
        auto p_earth = earth_opt.value().position.to_eigen_si();
        
        // 3. Manual Light-Time
        auto rho_geo = p_ast - p_earth;
        double tau_days = rho_geo.norm() / (constants::C_LIGHT * 1000.0 * 86400.0);
        auto t_emit = time::EpochTDB::from_mjd(t_obs.mjd() - tau_days);
        
        auto ast_emit_opt = horizons.query_vectors("234", t_emit, "@0");
        if (!ast_emit_opt) return 1;
        auto p_ast_emit = ast_emit_opt.value().position.to_eigen_si();
        
        auto rho_astro = p_ast_emit - p_earth;
        double r_astro = rho_astro.norm();
        double ra_calc = atan2(rho_astro(1), rho_astro(0)) * 180.0/M_PI;
        if (ra_calc < 0) ra_calc += 360.0;
        double dec_calc = asin(rho_astro(2) / r_astro) * 180.0/M_PI;
        
        cout << fixed << setprecision(10);
        cout << "\n--- RESULTS ---" << endl;
        cout << "JPL Quantity 1 RA:  " << obs_jpl.ra.to_deg() << " deg" << endl;
        cout << "JPL Quantity 1 Dec: " << obs_jpl.dec.to_deg() << " deg" << endl;
        cout << "Calculated RA:      " << ra_calc << " deg" << endl;
        cout << "Calculated Dec:     " << dec_calc << " deg" << endl;
        
        double dra = (ra_calc - obs_jpl.ra.to_deg()) * cos(obs_jpl.dec.to_rad()) * 3600.0;
        double ddec = (dec_calc - obs_jpl.dec.to_deg()) * 3600.0;
        
        cout << "\nDifference (Calc - JPL Q1):" << endl;
        cout << "dRA*cos(dec): " << dra * 1000.0 << " mas" << endl;
        cout << "dDec:         " << ddec * 1000.0 << " mas" << endl;

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
