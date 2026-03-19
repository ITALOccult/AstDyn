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
        
        // 1. Query JPL BARYCENTRIC VECTORS (CENTER='@0')
        cout << "Step: Querying JPL BARYCENTRIC VECTORS for (234) Barbara..." << endl;
        auto cart_opt = horizons.query_vectors("234", t_obs, "@0");
        if (!cart_opt) { cerr << "Failed to query vectors" << endl; return 1; }
        auto p_ast_ssb = cart_opt.value();

        // 2. Setup Astrometric Settings to match JPL Quantity 1
        astrometry::AstrometricSettings a_settings;
        a_settings.light_time_correction = true;
        a_settings.aberrazione_differenziale = false; 
        a_settings.deflessione_relativistica = false;  
        a_settings.frame_conversion_to_equatorial = true;

        // 3. Run Local Astrometry
        cout << "Step: Computing Astrometric Observation..." << endl;
        auto obs_res = astrometry::AstrometryReducer::compute_observation_from_cartesian(
            p_ast_ssb, t_obs, t_obs, cfg, a_settings);

        if (!obs_res) { cerr << "Astrometry failed" << endl; return 1; }
        auto obs = obs_res.value();

        cout << fixed << setprecision(10);
        cout << "ASTDYN Result (JPL Vectors + Local Astrometry):" << endl;
        cout << "RA:  " << obs.ra.to_hms() << " (" << obs.ra.to_deg() << " deg)" << endl;
        cout << "Dec: " << obs.dec.to_dms() << " (" << obs.dec.to_deg() << " deg)" << endl;

        // JPL Truth (Quantity 1 at 2026-May-04 00:00:00 UTC)
        double jpl_ra_deg = 197.9629235;
        double jpl_dec_deg = 13.7842222;

        double dra = (obs.ra.to_deg() - jpl_ra_deg) * cos(obs.dec.to_rad()) * 3600.0;
        double ddec = (obs.dec.to_deg() - jpl_dec_deg) * 3600.0;
        double sep = sqrt(dra*dra + ddec*ddec) * 1000.0;

        cout << "\nDifference (Local Astrometry - JPL Quantity 1):" << endl;
        cout << "dRA*cos(dec): " << dra << " arcsec (" << dra*1000.0 << " mas)" << endl;
        cout << "dDec:         " << ddec << " arcsec (" << ddec*1000.0 << " mas)" << endl;
        cout << "Total Sep:    " << sep << " mas" << endl;

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
