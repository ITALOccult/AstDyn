#ifndef ASTDYN_COORDINATES_FRAME_TRANSFORMS_HPP
#define ASTDYN_COORDINATES_FRAME_TRANSFORMS_HPP

#include "src/core/frame_tags.hpp"
#include "src/types/vectors.hpp"
#include "src/utils/time_types.hpp"
#include "rotation_matrices.hpp"

namespace astdyn::coordinates {

using core::GCRF;
using core::ITRF;
using types::Vector3;
using utils::Instant;

/** 
 * @brief Computes Greenwhich Mean Sidereal Time (GMST) for a given Instant.
 * Simplified for CTFYH example: linear approximation.
 */
[[nodiscard]] inline Radian compute_gmst(const Instant& epoch) noexcept {
    const double d_centuries = (epoch.mjd.value - 51544.5) / 36525.0;
    const double gmst_deg = 280.46061837 + 360.98564736629 * (epoch.mjd.value - 51544.5);
    return core::Radian(std::fmod(gmst_deg * (Constants::Math::PI / 180.0), Constants::Math::TWO_PI));
}

/** 
 * @brief Performs planetary rotation from GCRF (Celestial) to ITRF (Terrestrial).
 * Uses Sidereal Time approximation based on the epoch's Instant.
 */
template <typename Unit>
[[nodiscard]] constexpr Vector3<ITRF, Unit> rotate_gcrf_to_itrf(const Vector3<GCRF, Unit>& vec, const Instant& epoch) noexcept {
    const auto gmst = compute_gmst(epoch);
    const auto dcm = rotation_z(gmst);
    
    // Manual multiplication to preserve the new ITRF tag in return
    const double x = dcm.elements[0] * vec.x + dcm.elements[1] * vec.y;
    const double y = dcm.elements[3] * vec.x + dcm.elements[4] * vec.y;
    
    return Vector3<ITRF, Unit>(x, y, vec.z);
}

} // namespace astdyn::coordinates

#endif // ASTDYN_COORDINATES_FRAME_TRANSFORMS_HPP
