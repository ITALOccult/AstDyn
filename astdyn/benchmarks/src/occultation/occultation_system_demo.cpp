/**
 * @file occultation_system_demo.cpp
 * @brief Demo for searching occultations of a multi-body system (asteroid + satellites).
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/time/TimeScale.hpp"
#include <iostream>

using namespace astdyn;
using namespace astdyn::astrometry;

int main() {
    AstDynEngine engine;
    engine.set_verbose(true);

    // Window: March 2026
    time::EpochTDB start = time::EpochTDB::from_jd(time::calendar_to_mjd(2026, 3, 1, 0) + 2400000.5);
    time::EpochTDB end = time::EpochTDB::from_jd(time::calendar_to_mjd(2026, 3, 31, 23.9) + 2400000.5);

    // List of bodies in the system (NAIF IDs)
    // 87 = (87) Sylvia, 8701 = Romulus, 8702 = Remus (Common IDs for this system)
    std::vector<std::string> body_ids = {"87", "871", "872"}; 
    std::string bsp_path = "data/sylvia_system.bsp"; // Assuming this exists or the user will provide

    std::cout << "Searching for system occultations (Sylvia + Satellites)..." << std::endl;
    
    // In a real scenario, we'd use a real BSP. 
    // Here we wrap it in a try-catch because the file might be missing in this environment.
    try {
        auto system_occultations = OccultationLogic::find_system_occultations(
            body_ids, bsp_path, start, end, 12.0, engine
        );

        std::cout << "Found " << system_occultations.size() << " system occultation events." << std::endl;

        for (size_t i = 0; i < system_occultations.size(); ++i) {
            const auto& event = system_occultations[i];
            std::cout << "Event " << i << ": Star Gaia " << event.star.source_id 
                      << " Bodies involved: " << event.bodies.size() << std::endl;

            std::vector<OccultationPath> paths;
            std::vector<std::string> labels;
            for (const auto& body : event.bodies) {
                auto path = OccultationMapper::compute_path(
                    body.params,
                    event.star.ra, event.star.dec,
                    body.diameter,
                    time::to_utc(body.params.t_ca),
                    time::TimeDuration::from_seconds(600.0)
                );
                paths.push_back(path);
                labels.push_back(body.name);
            }

            std::string kml_name = "sylvia_system_event_" + std::to_string(i) + ".kml";
            OccultationMapper::export_kml(paths, labels, kml_name);
            std::cout << "  Exported: " << kml_name << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error during system search: " << e.what() << std::endl;
        std::cerr << "(This is expected if 'data/sylvia_system.bsp' is missing)" << std::endl;
    }

    return 0;
}
