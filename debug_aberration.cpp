#include "astdyn/AstDyn.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.ephemeris_type = EphemerisType::DE441;
        
        static AstDynEngine engine(cfg);
        static io::HorizonsClient horizons;
        
        double jd_utc = 2461164.5;
        auto t_utc = time::EpochTDB::from_jd(jd_utc); 
        
        std::cout << "Querying Horizons for Barbara (234) Geocentric (CENTER=500)..." << std::endl;
        // Table quantities: 1 (Astrometric RA), 3 (Apparent RA), 20 (Geometric RA)
        // We'll just use query_vectors to get the state and assume we have the right JD.
        auto state_opt = horizons.query_vectors("234", t_utc, "500");
        
        if (state_opt) {
            auto cart = state_opt.value(); 
            Eigen::Vector3d rho_eq = cart.position.to_eigen_si();
            double r = rho_eq.norm();
            
            auto ra_dec = [](const Eigen::Vector3d& v) {
                double r = v.norm();
                double ra = std::atan2(v(1), v(0));
                if (ra < 0) ra += constants::TWO_PI;
                double dec = std::asin(v(2) / r);
                return std::make_pair(ra, dec);
            };
            
            auto [ra_geom, dec_geom] = ra_dec(rho_eq);
            
            std::cout << std::fixed << std::setprecision(8);
            std::cout << "JPL VECTORS result (should be geometric/apparent depending on JPL internal state):\n";
            std::cout << "RA: " << (ra_geom * constants::RAD_TO_DEG / 15.0) << " h\n";
            
            // Apply Stellar Aberration
            auto de441 = std::make_shared<ephemeris::DE441Provider>(cfg.ephemeris_file);
            Eigen::Vector3d v_earth = de441->getVelocity(ephemeris::CelestialBody::EARTH, t_utc).to_eigen_si();
            
            Eigen::Vector3d p = rho_eq / r;
            Eigen::Vector3d beta = v_earth / (constants::C_LIGHT * 1000.0);
            double p_dot_v = p.dot(beta);
            double v2 = beta.squaredNorm();
            double inv_gamma = std::sqrt(1.0 - v2);
            double denom = 1.0 + p_dot_v;
            Eigen::Vector3d p_prime = (inv_gamma * p + (1.0 + p_dot_v / (1.0 + inv_gamma)) * beta) / denom;
            
            auto [ra_aber, dec_aber] = ra_dec(p_prime);
            std::cout << "RA with Stellar Aberration: " << (ra_aber * constants::RAD_TO_DEG / 15.0) << " h\n";
            
            std::cout << "\nJPL QUANTITIES (from my manual check):\n";
            std::cout << "Quantity 1 (Astrometric): 13.19753056 h (13h 11m 51.11s)\n";
            std::cout << "Quantity 3 (Apparent):    13.20163056 h (13h 12m 05.87s)\n";
            std::cout << "Quantity 20 (Geometric):  13.19782500 h (13h 11m 52.17s)\n";
            
            std::cout << "\nAnalysis:\n";
            std::cout << "Diff (Aber - JPL Astrometric): " << (ra_aber * constants::RAD_TO_DEG / 15.0 - 13.19753056) * 3600.0 << " seconds\n";
            std::cout << "Diff (Geom - JPL Astrometric): " << (ra_geom * constants::RAD_TO_DEG / 15.0 - 13.19753056) * 3600.0 << " seconds\n";
        }
    } catch (const std::exception& e) { std::cerr << "Error: " << e.what() << "\n"; }
    return 0;
}
