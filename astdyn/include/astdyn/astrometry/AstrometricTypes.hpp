#ifndef ASTDYN_ASTROMETRY_TYPES_HPP
#define ASTDYN_ASTROMETRY_TYPES_HPP

#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/core/physics_types.hpp"
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
    RightAscension ra;
    Declination dec;
    physics::Distance distance;
    Angle residual_ra;
    Angle residual_dec;

    AstrometricObservation() : ra(), dec(), distance(physics::Distance::from_au(1.0)), residual_ra(), residual_dec() {}
    AstrometricObservation(RightAscension r, Declination d, physics::Distance dist)
        : ra(r), dec(d), distance(dist), residual_ra(), residual_dec() {}
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
