#include "astdyn/astrometry/SPKChebyshevEphemeris.hpp"
#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/core/Constants.hpp"
#include <stdexcept>
#include <cmath>

namespace astdyn::astrometry {

SPKChebyshevEphemeris::SPKChebyshevEphemeris(
    int naif_id,
    io::SPKReader& reader,
    time::EpochTDB start_time,
    time::EpochTDB end_time,
    int degree)
    : start_epoch_(start_time), end_epoch_(end_time)
{
    jd_start_ = start_time.jd();
    jd_end_   = end_time.jd();

    // Create daily segments
    double current_jd = std::floor(jd_start_ - 0.5) + 0.5;
    while (current_jd < jd_end_) {
        time::EpochTDB midnight = time::EpochTDB::from_jd(current_jd);
        auto seg = catalog::fit_chebyshev_spk(reader, naif_id, midnight, 1.0, degree);
        segments_.push_back(seg);
        current_jd += 1.0;
    }
}

std::tuple<double, double, double> SPKChebyshevEphemeris::evaluate(time::EpochTDB epoch) const {
    return get_segment(epoch).evaluate_all(epoch.jd());
}

std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>> 
SPKChebyshevEphemeris::evaluate_full(time::EpochTDB epoch) const {
    return get_segment(epoch).evaluate_full(epoch.jd());
}

const catalog::ChebyshevSegment& SPKChebyshevEphemeris::get_segment(time::EpochTDB epoch) const {
    double jd = epoch.jd();
    if (jd < jd_start_ || jd > jd_end_) {
        throw std::out_of_range("SPKChebyshevEphemeris: epoch out of range");
    }

    size_t idx = static_cast<size_t>(std::floor(jd - std::floor(jd_start_ - 0.5) - 0.5));
    if (idx >= segments_.size()) {
        idx = segments_.size() - 1;
    }
    return segments_[idx];
}

} // namespace astdyn::astrometry
