#include "astdyn_wrapper.h"
#include <iostream>
#include <cmath>

namespace ioccultcalc {

/**
 * @brief Valida elementi Kepleriani prima della propagazione
 * @return true se validi, false altrimenti
 */
bool validateKeplerianElements(const astdyn::propagation::KeplerianElements& kep, const std::string& context = "") {
    bool valid = true;
    
    if (!context.empty()) {
        std::cout << "[VALIDATION] Context: " << context << std::endl;
    }
    
    // Check semi-major axis
    if (kep.semi_major_axis <= 0 || !std::isfinite(kep.semi_major_axis)) {
        std::cerr << "[VALIDATION] ❌ Invalid semi_major_axis: " << kep.semi_major_axis << std::endl;
        valid = false;
    }
    
    // Check eccentricity
    if (kep.eccentricity < 0 || kep.eccentricity >= 1.0 || !std::isfinite(kep.eccentricity)) {
        std::cerr << "[VALIDATION] ❌ Invalid eccentricity: " << kep.eccentricity << std::endl;
        valid = false;
    }
    
    // Check inclination
    if (kep.inclination < 0 || kep.inclination > M_PI || !std::isfinite(kep.inclination)) {
        std::cerr << "[VALIDATION] ❌ Invalid inclination: " << kep.inclination << " rad (" 
                  << kep.inclination * 180.0 / M_PI << " deg)" << std::endl;
        valid = false;
    }
    
    // Check longitude of ascending node
    if (!std::isfinite(kep.longitude_ascending_node)) {
        std::cerr << "[VALIDATION] ❌ Invalid longitude_ascending_node: " << kep.longitude_ascending_node << std::endl;
        valid = false;
    }
    
    // Check argument of perihelion
    if (!std::isfinite(kep.argument_perihelion)) {
        std::cerr << "[VALIDATION] ❌ Invalid argument_perihelion: " << kep.argument_perihelion << std::endl;
        valid = false;
    }
    
    // Check mean anomaly
    if (!std::isfinite(kep.mean_anomaly)) {
        std::cerr << "[VALIDATION] ❌ Invalid mean_anomaly: " << kep.mean_anomaly << std::endl;
        valid = false;
    }
    
    // Check epoch
    if (!std::isfinite(kep.epoch_mjd_tdb) || kep.epoch_mjd_tdb < 0) {
        std::cerr << "[VALIDATION] ❌ Invalid epoch_mjd_tdb: " << kep.epoch_mjd_tdb << std::endl;
        valid = false;
    }
    
    if (valid) {
        std::cout << "[VALIDATION] ✅ All elements valid" << std::endl;
    }
    
    return valid;
}

} // namespace ioccultcalc
