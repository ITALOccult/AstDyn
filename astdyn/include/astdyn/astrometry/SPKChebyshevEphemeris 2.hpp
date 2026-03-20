#ifndef ASTDYN_ASTROMETRY_SPK_CHEBYSHEV_EPHEMERIS_HPP
#define ASTDYN_ASTROMETRY_SPK_CHEBYSHEV_EPHEMERIS_HPP

#include "astdyn/astrometry/IChebyshevEphemeris.hpp"
#include "astdyn/io/SPKReader.hpp"
#include "astdyn/catalog/CatalogTypes.hpp"
#include <vector>
#include <memory>

namespace astdyn::astrometry {

class SPKChebyshevEphemeris : public IChebyshevEphemeris {
public:
    SPKChebyshevEphemeris(
        int naif_id,
        io::SPKReader& reader,
        time::EpochTDB start_time,
        time::EpochTDB end_time,
        int degree = 12
    );

    std::tuple<double, double, double> evaluate(time::EpochTDB epoch) const override;

    std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>> 
    evaluate_full(time::EpochTDB epoch) const override;

    time::EpochTDB start_epoch() const override { return start_epoch_; }
    time::EpochTDB end_epoch() const override { return end_epoch_; }
    size_t num_segments() const override { return segments_.size(); }
    const catalog::ChebyshevSegment& get_segment(time::EpochTDB epoch) const override;

private:
    std::vector<catalog::ChebyshevSegment> segments_;
    time::EpochTDB start_epoch_;
    time::EpochTDB end_epoch_;
    double jd_start_;
    double jd_end_;
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_SPK_CHEBYSHEV_EPHEMERIS_HPP
