#include "astdyn/AstDyn.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/astrometry/Astrometry.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace std;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        
        AstDynEngine engine;
        engine.set_config(cfg);
        
        // Star
        double star_ra_deg = 196.84357; 
        double star_dec_deg = 14.017908;
        
        // Orbit XML
        physics::KeplerianStateTyped<core::ECLIPJ2000> xml_orbit;
        xml_orbit.epoch = time::EpochTDB::from_jd(2461170.5);
        xml_orbit.a = physics::Distance::from_au(2.38449);
        xml_orbit.e = 0.24615;
        xml_orbit.i = astrometry::Angle::from_deg(15.3844);
        xml_orbit.node = astrometry::Angle::from_deg(144.4338);
        xml_orbit.omega = astrometry::Angle::from_deg(192.2935);
        xml_orbit.M = astrometry::Angle::from_deg(252.0926);
        xml_orbit.gm = physics::GravitationalParameter::sun();
        engine.set_initial_orbit(xml_orbit);

        double jd_val = 2461171.2078567; // The value from XML
        
        auto check_miss = [&](double jd, const string& label) {
            time::EpochTDB t = time::EpochTDB::from_jd(jd);
            auto pos = engine.propagate_to(t);
            astrometry::AstrometricSettings a_cfg;
            a_cfg.light_time_correction = true;
            a_cfg.aberrazione_differenziale = false; // Astrometric
            auto obs = *astrometry::AstrometryReducer::compute_observation(pos, xml_orbit.epoch, t, cfg, a_cfg);
            
            double d_ra = (obs.ra.to_deg() - star_ra_deg) * constants::DEG_TO_RAD;
            double d_dec = (obs.dec.to_deg() - star_dec_deg) * constants::DEG_TO_RAD;
            double xi = obs.distance.to_au() * cos(obs.dec.to_rad()) * sin(d_ra);
            double eta = obs.distance.to_au() * (sin(obs.dec.to_rad()) * cos(star_dec_deg * constants::DEG_TO_RAD) - cos(obs.dec.to_rad()) * sin(star_dec_deg * constants::DEG_TO_RAD) * cos(d_ra));
            double earth_rad_au = 6371.0 / 149597870.7;
            double dist = sqrt(xi*xi + eta*eta) / earth_rad_au;
            cout << label << " Dist: " << dist << " ER (" << (dist-1.0)*6371.0 << " km)" << endl;
        };

        cout << fixed << setprecision(10);
        check_miss(jd_val, "Assuming TDB=JD_VAL");
        check_miss(jd_val + 69.184/86400.0, "Assuming UTC=JD_VAL (TDB=JD_VAL+69s)");

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
