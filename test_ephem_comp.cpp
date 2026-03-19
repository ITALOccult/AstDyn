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
        
        // 1. Local Earth Position (DE441)
        auto de441 = make_shared<ephemeris::DE441Provider>(cfg.ephemeris_file);
        auto p_earth_local = de441->getState(ephemeris::CelestialBody::EARTH, t_obs).position.to_eigen_si();
        
        // 2. JPL Earth Position
        cout << "Querying JPL Earth Barycentric Vector..." << endl;
        auto earth_opt = horizons.query_vectors("399", t_obs, "@0");
        if (!earth_opt) return 1;
        auto p_earth_jpl = earth_opt.value().position.to_eigen_si();
        
        cout << fixed << setprecision(10);
        cout << "Earth Local (Bary, m): " << p_earth_local.transpose() << endl;
        cout << "Earth JPL   (Bary, m): " << p_earth_jpl.transpose() << endl;
        
        double diff = (p_earth_local - p_earth_jpl).norm();
        cout << "Earth Position Diff: " << diff << " meters" << endl;
        cout << "Ang. Error at 1.7 AU: " << (diff / (1.7 * 1.496e11)) * 206265.0 * 1000.0 << " mas" << endl;

        // 3. Check Sun position too
        auto p_sun_local = de441->getState(ephemeris::CelestialBody::SUN, t_obs).position.to_eigen_si();
        auto sun_opt = horizons.query_vectors("10", t_obs, "@0");
        if (!sun_opt) return 1;
        auto p_sun_jpl = sun_opt.value().position.to_eigen_si();
        cout << "Sun Position Diff: " << (p_sun_local - p_sun_jpl).norm() << " meters" << endl;

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
