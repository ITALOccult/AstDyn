/**
 * @file AsteroidChebyshevEphemeris.cpp
 * @brief Implementation of AsteroidChebyshevEphemeris.
 */

#include "astdyn/astrometry/AsteroidChebyshevEphemeris.hpp"
#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/AstDynEngine.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace astdyn::astrometry {

AsteroidChebyshevEphemeris::AsteroidChebyshevEphemeris(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
    time::EpochTDB start_time,
    time::EpochTDB end_time,
    const astdyn::AstDynConfig& config,
    int degree)
    : start_epoch_(start_time), end_epoch_(end_time)
{
    jd_start_ = start_time.jd();
    jd_end_ = end_time.jd();

    if (jd_end_ <= jd_start_) {
        throw std::invalid_argument("AsteroidChebyshevEphemeris: end_time must be after start_time.");
    }

    // Partition into daily segments (1 day = 1.0 JD)
    // We aim for segments starting at floor(jd) + 0.5 (noon) or floor(jd) (midnight)
    // For simplicity, we just divide the interval into 1-day chunks starting from jd_start_
    double current_jd = jd_start_;
    while (current_jd < jd_end_) {
        double next_jd = std::min(current_jd + 1.0, jd_end_);
        double segment_duration = next_jd - current_jd;
        
        // center_epoch for fit_chebyshev
        time::EpochTDB center = time::EpochTDB::from_jd(current_jd + segment_duration / 2.0);
        
        // Use the existing catalog utility to fit the segment
        auto seg = catalog::fit_chebyshev(initial_elements, center, segment_duration, config, degree);
        segments_.push_back(std::move(seg));
        
        current_jd = next_jd;
    }
}

std::tuple<double, double, double> AsteroidChebyshevEphemeris::evaluate(time::EpochTDB epoch) const {
    double jd = epoch.jd();
    
    if (jd < jd_start_ || jd > jd_end_) {
        throw std::out_of_range("AsteroidChebyshevEphemeris: requested epoch is outside the calculated interval.");
    }

    // Find the correct segment
    // Since segments are 1 day long starting from jd_start_:
    size_t index = static_cast<size_t>(std::floor(jd - jd_start_));
    if (index >= segments_.size()) {
        index = segments_.size() - 1; // Safeguard for jd == jd_end_
    }

    return segments_[index].evaluate_all(jd);
}

std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>> 
AsteroidChebyshevEphemeris::evaluate_full(time::EpochTDB epoch) const {
    double jd = epoch.jd();
    
    if (jd < jd_start_ || jd > jd_end_) {
        throw std::out_of_range("AsteroidChebyshevEphemeris: requested epoch is outside the calculated interval.");
    }

    size_t index = static_cast<size_t>(std::floor(jd - jd_start_));
    if (index >= segments_.size()) {
        index = segments_.size() - 1;
    }

    return segments_[index].evaluate_full(jd);
}

} // namespace astdyn::astrometry
