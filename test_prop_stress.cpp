#include "astdyn/AstDyn.hpp"
#include "astdyn/propagation/AASIntegrator.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace std;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.ephemeris_type = EphemerisType::DE441;
        cfg.integrator_type = IntegratorType::RKF78;
        cfg.tolerance = 1e-13;
        cfg.propagator_settings.include_asteroids = true;
        cfg.propagator_settings.use_default_asteroid_set = true;
        
        io::HorizonsClient horizons;
        auto t_start = time::EpochTDB::from_jd(2461164.5 - 1.0/24.0); // 1 hour before
        auto t_end = time::EpochTDB::from_jd(2461164.5);   // May 4, 2026 00:00
        
        // 1. Get JPL Truth at Start
        cout << "Fetching JPL State at 1 hour before target..." << endl;
        auto state_start_opt = horizons.query_vectors("234", t_start, "@0"); // Barycentric
        if (!state_start_opt) return 1;
        auto state_start = state_start_opt.value();

        // 2. Propagate with AstDyn
        cout << "Propagating 1 hour with RKF78 (1e-13) + Asteroids..." << endl;
        auto de441 = make_shared<ephemeris::DE441Provider>(cfg.ephemeris_file);
        auto ephem = make_shared<ephemeris::PlanetaryEphemeris>();
        ephem->setProvider(de441);
        
        auto integrator = make_shared<propagation::RKF78Integrator>(0.1, cfg.tolerance);
        propagation::Propagator prop(integrator, ephem, cfg.propagator_settings);
        
        auto state_end_calc = prop.propagate_cartesian(state_start, t_end);

        // 3. Get JPL Truth at End
        cout << "Fetching JPL State at May 4..." << endl;
        auto state_end_jpl_opt = horizons.query_vectors("234", t_end, "@0");
        if (!state_end_jpl_opt) return 1;
        auto state_end_jpl = state_end_jpl_opt.value();
        
        // 4. Compare
        auto p_calc = state_end_calc.position.to_eigen_si();
        auto p_jpl = state_end_jpl.position.to_eigen_si();
        double dist_err = (p_calc - p_jpl).norm();
        
        cout << fixed << setprecision(10);
        cout << "\n--- PROPAGATION RESULTS (3 DAYS) ---" << endl;
        cout << "Position Error: " << dist_err << " meters" << endl;
        
        // Convert to angular error as seen from Earth (~1.7 AU)
        double ang_err_mas = (dist_err / (1.7 * 1.496e11)) * 206265.0 * 1000.0;
        cout << "Angular Error (at 1.7 AU): " << ang_err_mas << " mas" << endl;

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
