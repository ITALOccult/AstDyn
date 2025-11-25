/**
 * @file Observation.cpp
 * @brief Implementation of observation utilities
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 */

#include "astdyn/observations/Observation.hpp"
#include <algorithm>
#include <stdexcept>

namespace astdyn {
namespace observations {

void ObservationSet::sortByTime() {
    std::sort(observations.begin(), observations.end(),
        [](const Observation& a, const Observation& b) {
            return a.getMJD() < b.getMJD();
        });
}

double ObservationSet::getTimeSpan() const {
    if (observations.empty()) return 0.0;
    
    double min_mjd = observations.front().getMJD();
    double max_mjd = observations.front().getMJD();
    
    for (const auto& obs : observations) {
        double mjd = obs.getMJD();
        if (mjd < min_mjd) min_mjd = mjd;
        if (mjd > max_mjd) max_mjd = mjd;
    }
    
    return max_mjd - min_mjd;
}

std::vector<OpticalObservation> ObservationSet::getOpticalObservations() const {
    std::vector<OpticalObservation> result;
    for (const auto& obs : observations) {
        if (obs.optical) {
            result.push_back(*obs.optical);
        }
    }
    return result;
}

std::vector<RadarObservation> ObservationSet::getRadarObservations() const {
    std::vector<RadarObservation> result;
    for (const auto& obs : observations) {
        if (obs.radar) {
            result.push_back(*obs.radar);
        }
    }
    return result;
}

CatalogCode parseCatalogCode(char code) {
    switch (code) {
        case 'a': return CatalogCode::USNO_A2;
        case 'b': return CatalogCode::USNO_B1;
        case 'c': return CatalogCode::UCAC1;
        case 'd': return CatalogCode::UCAC2;
        case 'e': return CatalogCode::UCAC3;
        case 'f': return CatalogCode::UCAC4;
        case 'g': return CatalogCode::USNO_SA2;
        case 'h': return CatalogCode::TYCHO_1;
        case 'i': return CatalogCode::TYCHO_2;
        case 'j': return CatalogCode::GSC_1_2;
        case 'k': return CatalogCode::GSC_2_3;
        case 'l': return CatalogCode::CMC_14;
        case 'm': return CatalogCode::CMC_15;
        case 'n': return CatalogCode::PPMXL;
        case 'o': return CatalogCode::GAIA_DR1;
        case 'p': return CatalogCode::GAIA_DR2;
        case 'q': return CatalogCode::GAIA_DR3;
        default: return CatalogCode::UNKNOWN;
    }
}

std::string catalogCodeToString(CatalogCode code) {
    switch (code) {
        case CatalogCode::USNO_A2: return "USNO-A2.0";
        case CatalogCode::USNO_B1: return "USNO-B1.0";
        case CatalogCode::UCAC1: return "UCAC-1";
        case CatalogCode::UCAC2: return "UCAC-2";
        case CatalogCode::UCAC3: return "UCAC-3";
        case CatalogCode::UCAC4: return "UCAC-4";
        case CatalogCode::USNO_SA2: return "USNO-SA2.0";
        case CatalogCode::TYCHO_1: return "Tycho-1";
        case CatalogCode::TYCHO_2: return "Tycho-2";
        case CatalogCode::GAIA_DR1: return "Gaia-DR1";
        case CatalogCode::GAIA_DR2: return "Gaia-DR2";
        case CatalogCode::GAIA_DR3: return "Gaia-DR3";
        case CatalogCode::GSC_1_2: return "GSC-1.2";
        case CatalogCode::GSC_2_3: return "GSC-2.3";
        case CatalogCode::CMC_14: return "CMC-14";
        case CatalogCode::CMC_15: return "CMC-15";
        case CatalogCode::PPMXL: return "PPMXL";
        default: return "Unknown";
    }
}

std::string observationTypeToString(ObservationType type) {
    switch (type) {
        case ObservationType::OPTICAL_RA_DEC: return "Optical RA/Dec";
        case ObservationType::OPTICAL_AZ_EL: return "Optical Az/El";
        case ObservationType::RADAR_RANGE: return "Radar Range";
        case ObservationType::RADAR_DOPPLER: return "Radar Doppler";
        case ObservationType::RADAR_RANGE_RATE: return "Radar Range-Rate";
        case ObservationType::PHOTOMETRY: return "Photometry";
        case ObservationType::OCCULTATION: return "Occultation";
        case ObservationType::SATELLITE_IMAGING: return "Satellite";
        case ObservationType::ROVING_OBSERVER: return "Roving";
        default: return "Unknown";
    }
}

} // namespace observations
} // namespace astdyn
