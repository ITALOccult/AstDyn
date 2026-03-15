#ifndef ASTDYN_ASTROMETRY_API_HPP
#define ASTDYN_ASTROMETRY_API_HPP

#include "astdyn/time/epoch.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "AstrometricTypes.hpp"
#include <Eigen/Dense>
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
     * @param initial Initial orbital state (Keplerian/Ecliptic).
     * @param t_elements Epoch of elements (time::EpochTDB).
     * @param t_obs Observation time (time::EpochTDB).
     * @param engine_cfg Configuration for the propagator.
     * @param a_cfg Settings for astrometric corrections.
     * @return AstrometricObservation or error.
     */
    static std::expected<AstrometricObservation, AstrometryError> compute_observation(
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial,
        const time::EpochTDB& t_elements,
        const time::EpochTDB& t_obs,
        const AstDynConfig& engine_cfg,
        const AstrometricSettings& a_cfg);

    /** @brief Computes an observation starting from truth Cartesian vector (meters, m/s). */
    static std::expected<AstrometricObservation, AstrometryError> compute_observation_from_cartesian(
        const physics::CartesianStateTyped<core::GCRF>& initial,
        const time::EpochTDB& t_elements,
        const time::EpochTDB& t_obs,
        const AstDynConfig& engine_cfg,
        const AstrometricSettings& a_cfg);

    static Eigen::Vector3d apply_stellar_aberration(
        const Eigen::Vector3d& rho_eq,
        const Eigen::Vector3d& earth_velocity_eq);

    /**
     * @brief Compute the 1-sigma cross-track uncertainty in the B-plane.
     * 
     * @param covariance_eq 6x6 covariance matrix in Equatorial J2000 (km, km/s)
     * @param rho_eq Vector from observer to body (meters)
     * @param velocity_eq Velocity of body (m/s)
     * @param star_ra Star RA
     * @param star_dec Star Dec
     * @return 1-sigma cross-track uncertainty in km
     */
    static double compute_cross_track_uncertainty(
        const Eigen::Matrix<double, 6, 6>& covariance_eq,
        const Eigen::Vector3d& rho_eq,
        const Eigen::Vector3d& velocity_eq,
        double star_ra_rad,
        double star_dec_rad);

    /** @brief Step 2b: Gravitational Light Deflection (Sun). */
    static Eigen::Vector3d apply_light_deflection(
        const Eigen::Vector3d& rho_eq,
        const Eigen::Vector3d& earth_to_sun_eq);

private:
    /** @brief Step 1: Iterative Light-Time Correction Kernel. */
    static Eigen::Vector3d compute_light_time_corrected_pos(
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial,
        const time::EpochTDB& t_elements,
        const time::EpochTDB& t_obs,
        const Eigen::Vector3d& earth_pos_helio_ecl,
        const AstDynConfig& cfg);

    /** @brief Step 3: Frame Transformation (Ecliptic -> Equatorial). */
    static Eigen::Vector3d convert_frame_if_needed(
        const Eigen::Vector3d& vec,
        const AstrometricSettings& settings);

    /** @brief Step 4: Final RA/Dec/Distance conversion. */
    static AstrometricObservation finalize_observation(
        const Eigen::Vector3d& final_rho_eq);
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_API_HPP
