#ifndef ASTDYN_COORDINATES_FRAME_TRANSFORMS_HPP
#define ASTDYN_COORDINATES_FRAME_TRANSFORMS_HPP

#include "astdyn/core/frame_tags.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "astdyn/time/epoch.hpp"
#include "rotation_matrices.hpp"

namespace astdyn::coordinates {

using core::GCRF;
using core::ITRF;
using math::Vector3;

/** 
 * @brief Computes Greenwhich Mean Sidereal Time (GMST) for a given Epoch.
 * Simplified for CTFYH example: linear approximation.
 */
[[nodiscard]] inline Radian compute_gmst(const time::EpochUTC& epoch) noexcept {
    const double d_centuries = (epoch.mjd() - 51544.5) / 36525.0;
    const double gmst_deg = 280.46061837 + 360.98564736629 * (epoch.mjd() - 51544.5);
    return core::Radian(std::fmod(gmst_deg * (Constants::Math::PI / 180.0), Constants::Math::TWO_PI));
}

/** 
 * @brief Performs planetary rotation from GCRF (Celestial) to ITRF (Terrestrial).
 * Uses Sidereal Time approximation based on the epoch's time.
 */
template <typename Unit>
[[nodiscard]] constexpr Vector3<ITRF, Unit> rotate_gcrf_to_itrf(const Vector3<GCRF, Unit>& vec, const time::EpochUTC& epoch) noexcept {
    const auto gmst = compute_gmst(epoch);
    const auto dcm = rotation_z(gmst);
    
    // Manual multiplication to preserve the new ITRF tag in return
    const double x = dcm.elements[0] * vec.x_si() + dcm.elements[1] * vec.y_si();
    const double y = dcm.elements[3] * vec.x_si() + dcm.elements[4] * vec.y_si();
    
    return Vector3<ITRF, Unit>::from_si(x, y, vec.z_si());
}

} // namespace astdyn::coordinates

#endif // ASTDYN_COORDINATES_FRAME_TRANSFORMS_HPP
