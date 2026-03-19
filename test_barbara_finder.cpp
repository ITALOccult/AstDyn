#include "astdyn/AstDyn.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/astrometry/ClosestApproachFinder.hpp"
#include "astdyn/astrometry/ChebyshevEphemerisManager.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::astrometry;
using namespace std;

int main() {
    try {
        // 1. Setup Engine
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.integrator_type = IntegratorType::RKF78;
        cfg.verbose = false;
        
        AstDynEngine engine;
        engine.set_config(cfg);
        
        // 2. Initialize Gaia Catalog
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa"})");
        
        // 3. Get Star (Gaia DR3 3737957664602055808)
        auto star_opt = catalog::GaiaDR3Catalog::instance().by_source_id(3737957664602055808);
        if (!star_opt) {
            cerr << "Error: Star not found in Gaia catalog." << endl;
            return 1;
        }
        const auto& star = *star_opt;
        cout << "Target Star: " << star.source_id << " (G=" << star.g_mag << ")" << endl;
        cout << "Star position (ICRS): RA=" << star.ra.to_deg() << " deg, Dec=" << star.dec.to_deg() << " deg" << endl;

        // 4. Get JPL Nominal Orbit for Barbara (234)
        io::HorizonsClient horizons;
        string id = "234";
        time::EpochTDB t_start = time::EpochTDB::from_jd(2461171.0); // May 10, 2026 12:00 UTC
        time::EpochTDB t_end   = time::EpochTDB::from_jd(2461171.4); // May 10, 2026 21:30 UTC
        
        cout << "Querying JPL Horizons for Barbara (234) nominal elements..." << endl;
        auto elements = horizons.query_elements(id, t_start);
        if (!elements) {
            cerr << "Error: Failed to fetch JPL elements." << endl;
            return 1;
        }

        auto physical = horizons.query_physical_properties(id);
        double diameter = physical ? physical->diameter_km : 45.9;
        cout << "Asteroid Diameter: " << diameter << " km" << endl;

        // 5. Generate Chebyshev Ephemeris
        cout << "Generating Chebyshev ephemeris..." << endl;
        ChebyshevEphemerisManager manager(cfg);
        manager.add_asteroid(id, *elements, t_start, t_end);
        const auto& ephem = manager.get_ephemeris(id);

        // 6. Find Closest Approach using the NEW Finder
        cout << "Running ClosestApproachFinder..." << endl;
        // Define test_jd and config for the new API call
        double test_jd = t_start.jd(); // Or any JD within the ephemeris segment
        ClosestApproachFinder::Config config; // Default config
        config.max_shadow_distance = Angle::from_arcsec(120.0); // Example value

        auto approaches = ClosestApproachFinder::find_in_segment(
            ephem.get_segment(test_jd), star, t_start, t_end, 
            config.max_shadow_distance);

        cout << fixed << setprecision(6);
        cout << "Found " << approaches.size() << " results." << endl;
        for (const auto& ca : approaches) {
            cout << "------------------------------------------------" << endl;
            cout << "T_CA (TDB): " << ca.t_ca.jd() << endl;
            // Calculate asteroid RA/Dec at TCA
            auto [p_tca, v_tca] = ephem.get_segment(ca.t_ca).evaluate_full(ca.t_ca.jd());
            cout << "Asteroid at TCA: RA=" << std::get<0>(p_tca) << " deg, Dec=" << std::get<1>(p_tca) << " deg" << endl;
            cout << "Separation: " << ca.separation.to_arcsec() << " arcsec" << endl;
            cout << "PA:         " << ca.position_angle_deg << " deg" << endl;
            cout << "Rate:       " << ca.relative_rate_arcsec_min << " arcsec/min" << endl;
            cout << "App Diam:   " << ca.apparent_diameter_as << " arcsec" << endl;
            
            bool is_hit = ca.impact_parameter < core::Distance::from_km(6371.0 + diameter/2.0);
            cout << "Occultation: " << (is_hit ? "YES" : "NO") << endl;
            cout << "Impact (Geocenter): " << ca.impact_parameter.to_km() << " km"
                 << " (" << ca.impact_parameter.to_km() / 6371.0 << " Earth Radii)" << endl;
            if (is_hit) cout << "=> SHADOW HITS EARTH!" << endl;
            else cout << "=> SHADOW MISSES EARTH (GEOCENTRIC)" << endl;
        }

        // 7. Test with User XML Elements
        cout << "\n==============================================" << endl;
        cout << "Testing with USER XML ELEMENTS..." << endl;
        
        physics::KeplerianStateTyped<core::ECLIPJ2000> xml_elements;
        xml_elements.a = physics::Distance::from_au(2.38449);
        xml_elements.e = 0.24615;
        xml_elements.i = Angle::from_deg(15.3844);
        xml_elements.node = Angle::from_deg(144.4338);
        xml_elements.omega = Angle::from_deg(192.2935);
        xml_elements.M = Angle::from_deg(252.0926);
        xml_elements.epoch = time::EpochTDB::from_jd(2461170.5 + 0.2926); // Decided by M? No, usually t_0.
        // The XML says 2026,5,10. JD for 2026-05-10 00:00 is 2461170.5
        xml_elements.epoch = time::EpochTDB::from_jd(2461170.5);

        ChebyshevEphemerisManager xml_manager(cfg);
        xml_manager.add_asteroid(id + "_xml", xml_elements, t_start, t_end);
        const auto& xml_ephem = xml_manager.get_ephemeris(id + "_xml");

        auto xml_approaches = ClosestApproachFinder::find(
            xml_ephem, *star, t_start, t_end, 
            Angle::from_arcsec(120.0), 1440, diameter);

        cout << "Found " << xml_approaches.size() << " results." << endl;

        for (const auto& ca : xml_approaches) {
            cout << "------------------------------------------------" << endl;
            cout << "T_CA (TDB): " << ca.t_ca.jd() << endl;
            cout << "Separation: " << ca.separation.to_arcsec() << " arcsec" << endl;
            auto [pos, vel] = xml_ephem.get_segment(ca.t_ca).evaluate_full(ca.t_ca.jd());
            double dist_au = std::get<2>(pos);
            double impact_km = ca.separation.to_rad() * dist_au * 1.495978707e8;
            double er = impact_km / 6371.0;
            
            cout << "Impact (Geocenter): " << impact_km << " km (" << er << " Earth Radii)" << endl;
            if (er < 1.0) {
                cout << "=> SHADOW AXES HITS EARTH DISK (GEOCENTRIC)" << endl;
            } else {
                cout << "=> SHADOW MISSES EARTH (GEOCENTRIC)" << endl;
            }
        }

    } catch (const exception& e) {
        cerr << "Exception: " << e.what() << endl;
        return 1;
    }
    return 0;
}
