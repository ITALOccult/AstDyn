#include "astdyn/AstDyn.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;
using namespace std;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.ephemeris_type = EphemerisType::DE441;
        cfg.integrator_type = IntegratorType::RKF78; 
        cfg.tolerance = 1e-12;
        cfg.propagator_settings.include_planets = true;
        cfg.propagator_settings.include_asteroids = true;
        cfg.propagator_settings.use_default_asteroid_set = true;
        
        static AstDynEngine engine(cfg);
        
        double jd_utc = 2461164.5;
        auto t_utc = time::EpochTDB::from_jd(jd_utc); 
        
        io::HorizonsClient horizons;
        auto kep_opt = horizons.query_elements("234", t_utc); 
        
        if (!kep_opt) {
            cerr << "Failed to query elements from Horizons" << endl;
            return 1;
        }
        
        auto elements = kep_opt.value();
        
        // Match Quantity 1: No Stellar Aberration, No Light Deflection
        cfg.aberrazione_differenziale = false;
        cfg.deflessione_relativistica = false;
        
        auto obs_res = astrometry::AstrometryReducer::compute_topocentric_observation(
            elements, t_utc, t_utc, "500", cfg
        );
        
        if (obs_res) {
            auto obs = obs_res.value();
            cout << fixed << setprecision(10);
            cout << "ASTDYN Result (Aberration OFF, GCRF/ICRF):\n";
            cout << "RA:  " << obs.ra.to_hms() << " (" << obs.ra.to_deg() << " deg)\n";
            cout << "Dec: " << obs.dec.to_dms() << " (" << obs.dec.to_deg() << " deg)\n";
            
            // JPL EXPECTED (Quantity 1 at 2026-05-04 00:00:00 UT)
            // 2026-May-04 00:00     13 11 51.10 +13 47 03.2
            double jpl_ra_deg = (13.0 + 11.0/60.0 + 51.10/3600.0) * 15.0; 
            double jpl_dec_deg = 13.0 + 47.0/60.0 + 03.2/3600.0;
            
            double dra = (obs.ra.to_deg() - jpl_ra_deg) * 3600.0 * cos(obs.dec.to_rad());
            double ddec = (obs.dec.to_deg() - jpl_dec_deg) * 3600.0;
            
            cout << "\nDifference (AstDyn - JPL Quantity 1):\n";
            cout << "dRA*cos(dec): " << dra << " arcsec (" << dra*1000.0 << " mas)\n";
            cout << "dDec:         " << ddec << " arcsec (" << ddec*1000.0 << " mas)\n";
            cout << "Total Sep:    " << sqrt(dra*dra + ddec*ddec) * 1000.0 << " mas\n";
        } else {
            cerr << "Astrometry calculation failed" << endl;
        }
        
    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
