#include "astdyn/AstDyn.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace std;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.ephemeris_type = EphemerisType::DE441;
        
        io::HorizonsClient horizons;
        auto t_obs = time::EpochTDB::from_jd(2461164.5); // May 4, 2026
        
        // 1. Query JPL BARYCENTRIC VECTORS at t_obs
        cout << "Querying JPL BARYCENTRIC VECTORS at t_obs..." << endl;
        auto ast_opt = horizons.query_vectors("234", t_obs, "@0");
        auto earth_opt = horizons.query_vectors("399", t_obs, "@0");
        if (!ast_opt || !earth_opt) return 1;
        
        auto p_ast = ast_opt.value().position.to_eigen_si();
        auto p_earth = earth_opt.value().position.to_eigen_si();
        
        // 2. RAW GEOMETRIC RA/DEC (No light-time)
        auto rho_geo = p_ast - p_earth;
        double r_geo = rho_geo.norm();
        double ra_geo = atan2(rho_geo(1), rho_geo(0)) * 180.0/M_PI;
        if (ra_geo < 0) ra_geo += 360.0;
        double dec_geo = asin(rho_geo(2) / r_geo) * 180.0/M_PI;
        
        cout << fixed << setprecision(10);
        cout << "Raw Geometric (Local, No Light-Time):" << endl;
        cout << "RA:  " << ra_geo << " deg" << endl;
        cout << "Dec: " << dec_geo << " deg" << endl;

        // 3. Query JPL VECTORS at t_emit (t_obs - tau)
        double tau_days = r_geo / (constants::C_LIGHT * 1000.0 * 86400.0);
        auto t_emit = time::EpochTDB::from_mjd(t_obs.mjd() - tau_days);
        cout << "\nLight-time tau: " << tau_days * 86400.0 << " seconds (" << tau_days << " days)" << endl;
        
        auto ast_emit_opt = horizons.query_vectors("234", t_emit, "@0");
        if (!ast_emit_opt) return 1;
        auto p_ast_emit = ast_emit_opt.value().position.to_eigen_si();
        
        // 4. ASTROMETRIC RA/DEC (With light-time)
        auto rho_astro = p_ast_emit - p_earth;
        double r_astro = rho_astro.norm();
        double ra_astro = atan2(rho_astro(1), rho_astro(0)) * 180.0/M_PI;
        if (ra_astro < 0) ra_astro += 360.0;
        double dec_astro = asin(rho_astro(2) / r_astro) * 180.0/M_PI;
        
        cout << "Astrometric (Local, with JPL Emit Vector):" << endl;
        cout << "RA:  " << ra_astro << " deg" << endl;
        cout << "Dec: " << dec_astro << " deg" << endl;

        // 5. Compare with JPL Quantity 1
        double jpl_ra_deg = 197.9629235;
        double jpl_dec_deg = 13.7842222;
        cout << "\nJPL Quantity 1 (Astrometric):" << endl;
        cout << "RA:  " << jpl_ra_deg << " deg" << endl;
        cout << "Dec: " << jpl_dec_deg << " deg" << endl;

        double err_astro = (ra_astro - jpl_ra_deg) * cos(dec_astro * M_PI/180.0) * 3600.0;
        cout << "Error in RA*cos(dec): " << err_astro * 1000.0 << " mas" << endl;

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
