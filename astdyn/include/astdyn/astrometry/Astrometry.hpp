#ifndef ASTDYN_ASTROMETRY_API_HPP
#define ASTDYN_ASTROMETRY_API_HPP

#include "astdyn/time/epoch.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "AstrometricTypes.hpp"
#include <Eigen/Dense>
#include <expected>

#include <memory>

namespace astdyn {
    class AstDynEngine;
    namespace ephemeris { class DE441Provider; }
}

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

    /**
     * @brief Computes a topocentric RA/Dec observation from an observatory code.
     */
    static std::expected<AstrometricObservation, AstrometryError> compute_topocentric_observation(
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial,
        const time::EpochTDB& t_elements,
        const time::EpochTDB& t_obs,
        const std::string& obs_code,
        const AstDynConfig& engine_cfg);

    /**
     * @brief Applica l'aberrazione differenziale (Stellar/Planetary Aberration).
     * Deve essere applicata DOPO la correzione del tempo di luce.
     * Corregge lo spostamento apparente dovuto al moto dell'osservatore durante il tragitto dei fotoni.
     */
    static Eigen::Vector3d aberrazione_differenziale(
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

    /**
     * @brief Applica la deflessione relativistica (Gravitational Light Deflection).
     * Deve essere applicata su coordinate già corrette per il tempo di luce.
     * Modella la curvatura dello spazio-tempo causata dalla massa del Sole (effetto Schwarzschild).
     */
    static Eigen::Vector3d deflessione_relativistica(
        const Eigen::Vector3d& rho_eq,
        const Eigen::Vector3d& earth_to_sun_eq);

private:
    static std::shared_ptr<::astdyn::ephemeris::DE441Provider> sync_ephemeris(const std::string& path);

    static Eigen::Vector3d compute_earth_helio(std::shared_ptr<::astdyn::ephemeris::DE441Provider> de441, const time::EpochTDB& t_obs);

    static Eigen::Vector3d compute_light_time_corrected_pos(
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial,
        const time::EpochTDB& t_elements,
        const time::EpochTDB& t_obs,
        const Eigen::Vector3d& earth_pos_helio_ecl,
        const AstDynConfig& cfg);

    static AstrometricObservation finalize_observation(const Eigen::Vector3d& final_rho_eq);
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_API_HPP
