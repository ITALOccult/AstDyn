#include "astdyn/AstDyn.hpp"
#include "astdyn/astrometry/Astrometry.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;
using namespace std;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.integrator_type = IntegratorType::RKF78;
        
        astrometry::AstrometricSettings a_cfg;
        a_cfg.light_time_correction = true;
        a_cfg.aberrazione_differenziale = false; // JPL Astrometric doesn't include stellar aberration
        a_cfg.deflessione_relativistica = true;
        
        // Star coordinates (Gaia DR3 predicted at 2026.3)
        double star_ra = 197.924349;
        double star_dec = 13.832289;
        
        // Barbara orbit (May 4.0 00:00 UTC)
        time::EpochTDB t0 = time::EpochTDB::from_jd(2461164.5);
        io::HorizonsClient horizons;
        auto state_res = horizons.query_elements("234", t0);
        if (!state_res) return 1;
        auto initial = *state_res;
        
        cout << fixed << setprecision(6);
        cout << "Time (UTC), JD, Sep (arcsec), Shadow_Dist (km)" << endl;
        
        // Scan around 03:30 UTC
        for (int m = 0; m <= 60; m += 1) {
            double jd = 2461164.5 + (3.0/24.0) + (m / 1440.0);
            time::EpochTDB t = time::EpochTDB::from_jd(jd);
            
            auto obs_res = astrometry::AstrometryReducer::compute_observation(initial, t0, t, cfg, a_cfg);
            if (!obs_res) continue;
            auto obs = *obs_res;
            
            double d_ra = (obs.ra.to_deg() - star_ra) * cos(obs.dec.to_rad());
            double d_dec = obs.dec.to_deg() - star_dec;
            double sep_deg = sqrt(d_ra*d_ra + d_dec*d_dec);
            double sep_arcsec = sep_deg * 3600.0;
            
            // Distance from Earth center to shadow axis
            double dist_km = sep_deg * 0.0174532925 * obs.distance.to_km();
            
            cout << "03:" << (m < 10 ? "0" : "") << m << ":00, " 
                 << jd << ", " << sep_arcsec << ", " << dist_km << endl;
        }

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
