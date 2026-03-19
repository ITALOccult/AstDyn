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
        
        double jd_utc = 2461164.5;
        auto t_utc = time::EpochTDB::from_jd(jd_utc); 
        
        io::HorizonsClient horizons;
        
        // 1. Query GEOMETRIC VECTORS from JPL (Geocentric GCRF)
        auto vec_opt = horizons.query_vectors("234", t_utc, "500");
        if (!vec_opt) return 1;
        auto jpl_geoc_gcrf = vec_opt.value().position.to_eigen_si();
        
        // 2. Query HELIOCENTRIC ELEMENTS from JPL
        auto kep_opt = horizons.query_elements("234", t_utc);
        if (!kep_opt) return 1;
        auto elements = kep_opt.value();
        
        // 3. Convert elements to Cartesian Heliocentric
        auto cart_helio_ecl = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(elements);
        auto p_ast_helio_ecl = cart_helio_ecl.position.to_eigen_si();
        
        // Transform Heliocentric Ecliptic to Heliocentric GCRF
        auto p_ast_helio_gcrf = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(
            math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(p_ast_helio_ecl.x(), p_ast_helio_ecl.y(), p_ast_helio_ecl.z())
        ).to_eigen_si();
        
        // 4. Get Earth position from DE441 (Heliocentric)
        auto de441 = std::make_shared<ephemeris::DE441Provider>(cfg.ephemeris_file);
        auto p_earth_bary = de441->getPosition(ephemeris::CelestialBody::EARTH, t_utc).to_eigen_si();
        auto p_sun_bary = de441->getPosition(ephemeris::CelestialBody::SUN, t_utc).to_eigen_si();
        Eigen::Vector3d p_earth_helio_gcrf = p_earth_bary - p_sun_bary;
        
        // 5. Calculate Geocentric Position
        // Pos(Ast)_Geoc = Pos(Ast)_Helio - Pos(Earth)_Helio
        Eigen::Vector3d p_ast_geoc_gcrf = p_ast_helio_gcrf - p_earth_helio_gcrf;
        
        cout << fixed << setprecision(10);
        cout << "JPL Geocentric (GCRF): " << jpl_geoc_gcrf.transpose() << endl;
        cout << "Local Geocentric (GCRF): " << p_ast_geoc_gcrf.transpose() << endl;
        
        double diff = (jpl_geoc_gcrf - p_ast_geoc_gcrf).norm();
        cout << "Difference magnitude: " << diff / 1000.0 << " km" << endl;
        
        // Convert back to RA/Dec
        auto ra_dec = [](const Eigen::Vector3d& v) {
            double r = v.norm();
            double ra = std::atan2(v(1), v(0));
            if (ra < 0) ra += constants::TWO_PI;
            double dec = std::asin(v(2) / r);
            return std::make_pair(ra, dec);
        };
        auto [ra, dec] = ra_dec(p_ast_geoc_gcrf);
        cout << "Local RA:  " << (ra * constants::RAD_TO_DEG / 15.0) << " h" << endl;

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
