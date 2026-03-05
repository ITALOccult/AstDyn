#ifndef ASTDYN_ASTROMETRY_API_HPP
#define ASTDYN_ASTROMETRY_API_HPP

#include "src/utils/time_types.hpp"
#include "src/core/units.hpp"
#include "src/core/frame_tags.hpp"
#include "astdyn/propagation/Propagator.hpp"

#include "AstrometricTypes.hpp"
#include "src/types/orbital_state.hpp"
#include "astdyn/AstDynEngine.hpp"
#include <Eigen/Dense>
#include <optional>

namespace astdyn {
    struct AstDynConfig;
}

#include <expected>

namespace astdyn::astrometry {

/**
 * @brief Top-level API for computing high-precision astrometric observations.
 */
class AstrometryReducer {
public:
    /**
     * @brief Computes an astrometric RA/Dec observation with corrections.
     * 
     * @param initial_state Initial orbital state (Keplerian/Ecliptic).
     * @param t_obs Observation time (Instant).
     * @param p_cfg Configuration for the propagator.
     * @param a_cfg Settings for astrometric corrections.
     * @return AstrometricObservation or error.
     */
    static std::expected<AstrometricObservation, AstrometryError> compute_observation(
        const types::OrbitalState<core::ECLIPJ2000, types::KeplerianTag>& initial_state,
        const utils::Instant& t_obs,
        const AstDynConfig& engine_cfg,
        const AstrometricSettings& a_cfg);

private:
    /** @brief Step 1: Iterative Light-Time Correction Kernel. */
    static Eigen::Vector3d compute_light_time_corrected_pos(
        const types::OrbitalState<core::ECLIPJ2000, types::KeplerianTag>& initial,
        const utils::Instant& t_obs,
        const Eigen::Vector3d& earth_pos,
        const AstDynConfig& cfg);

    /** @brief Step 2: Annual Stellar Aberration Correction. */
    static Eigen::Vector3d apply_stellar_aberration(
        const Eigen::Vector3d& geometric_rho,
        const Eigen::Vector3d& earth_velocity);

    /** @brief Step 3: Frame Transformation (Ecliptic -> Equatorial). */
    static Eigen::Vector3d convert_frame_if_needed(
        const Eigen::Vector3d& vec,
        const AstrometricSettings& settings);

    /** @brief Step 4: Final RA/Dec/Distance conversion. */
    static AstrometricObservation finalize_observation(
        const Eigen::Vector3d& final_rho);
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_API_HPP
