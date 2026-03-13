#ifndef ASTDYN_ASTROMETRY_TYPES_HPP
#define ASTDYN_ASTROMETRY_TYPES_HPP

#include "src/core/units.hpp"
#include "astdyn/propagation/Propagator.hpp"

namespace astdyn::astrometry {

/**
 * @brief Configuration for high-precision astrometric reduction.
 */
struct AstrometricSettings {
    bool light_time_correction = true;
    bool stellar_aberration = true;
    bool light_deflection = true;
    bool frame_conversion_to_equatorial = true;
};

/**
 * @brief Result of an astrometric observation computation.
 */
struct AstrometricObservation {
    core::Radian ra;
    core::Radian dec;
    core::Meter distance;
    core::MilliArcSecond residual_ra{0.0};
    core::MilliArcSecond residual_dec{0.0};
};

/**
 * @brief Error codes for astrometry module.
 */
enum class AstrometryError {
    PropagationFailed,
    EphemerisUnavailable,
    ConvergenceFailed,
    InvalidInput
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_TYPES_HPP
