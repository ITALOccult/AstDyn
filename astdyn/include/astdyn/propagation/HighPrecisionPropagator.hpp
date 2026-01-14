/**
 * @file HighPrecisionPropagator.hpp
 * @brief Easy-to-use API for high-precision orbit propagation and geocentric observation
 * @author AstDyn Team
 */

#ifndef ASTDYN_HIGH_PRECISION_PROPAGATOR_HPP
#define ASTDYN_HIGH_PRECISION_PROPAGATOR_HPP

#include "astdyn/propagation/OrbitalElements.hpp"
#include <string>
#include <memory>
#include <Eigen/Dense>

namespace astdyn::ephemeris {
    class PlanetaryEphemeris;
    class EphemerisProvider;
}

namespace astdyn::propagation {

class Propagator;

/**
 * @brief High-level API for precision propagation
 * 
 * Encapsulates:
 * - High-fidelity N-body propagation settings
 * - JPL DE441 ephemeris loading (optional)
 * - Light-time correction
 * - Coordinate frame handling (Ecliptic <-> Equatorial)
 */
class HighPrecisionPropagator {
public:
    enum class InputFrame {
        ECLIPTIC,   ///< Mean Ecliptic J2000 (standard for AstDys/MPC)
        EQUATORIAL  ///< J2000 Equatorial / ICRF
    };

    struct Config {
        std::string de441_path = "";     ///< Path to de44xx.bsp. If empty, uses analytical ephemeris.
        std::string asteroid_ephemeris_file = ""; ///< Path to asteroid SPK kernel
        double step_size;     ///< Integrator initial step [days]
        double tolerance;   ///< Integrator tolerance
        
        // Perturbations
        bool perturbations_planets;
        bool perturbations_asteroids;
        bool relativity;

        Config() 
            : step_size(0.5), tolerance(1e-13), 
              perturbations_planets(true), perturbations_asteroids(true), relativity(true) {}
    };

    struct ObservationResult {
        double ra_deg;              ///< Right Ascension J2000 [deg]
        double dec_deg;             ///< Declination J2000 [deg]
        double distance_au;         ///< Geometric distance [AU]
        double light_time_sec;      ///< Light travel time [s]
        Eigen::Vector3d geocentric_position; ///< Vector Earth->Body (J2000 Equatorial) [AU]
    };

    /**
     * @brief Initialize propagator with configuration
     * 
     * @example
     * HighPrecisionPropagator::Config config;
     * config.de441_path = "path/to/de441.bsp";
     * config.perturbations_planets = true;
     * 
     * HighPrecisionPropagator propagator(config);
     * auto result = propagator.calculateGeocentricObservation(elements, target_jd);
     * 
     * @param config Configuration settings for the propagator
     */
    explicit HighPrecisionPropagator(const Config& config = Config());
    ~HighPrecisionPropagator();

    /**
     * @brief Calculate geocentric astrometric position
     * 
     * Performs:
     * 1. N-Body propagation from initial epoch to target time
     * 2. Light-time correction loop
     * 3. Precise Earth position retrieval
     * 4. Conversion to RA/DEC J2000
     * 
     * @param initial_elements Osculating elements at starting epoch (must be valid Mean or Osculating depending on usage, but usually Osculating J2000 Ecliptic)
     * @param target_jd_tdb Target time (JD TDB)
     * @return ObservationResult
     */
    ObservationResult calculateGeocentricObservation(
        const KeplerianElements& initial_elements, 
        double target_jd_tdb,
        InputFrame frame = InputFrame::ECLIPTIC
    );

    /**
     * @brief Get underlying planetary ephemeris provider
     */
    std::shared_ptr<astdyn::ephemeris::EphemerisProvider> getEphemerisProvider() const;

    /**
     * @brief Create a fresh configured Propagator instance
     */
    std::unique_ptr<Propagator> createPropagator() const;

    /**
     * @brief High-precision propagation to target epoch
     * 
     * @param initial_elements Initial orbit
     * @param target_mjd_tdb Target epoch
     * @param frame Input frame
     * @return Cartesian state at target epoch (ICRF Equatorial)
     */
    CartesianElements propagate_cartesian(
        const KeplerianElements& initial_elements,
        double target_mjd_tdb,
        InputFrame frame = InputFrame::ECLIPTIC
    );

private:
private:
    Config config_;
    
    // Shared resources
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> planetary_ephemeris_;
    std::shared_ptr<astdyn::ephemeris::EphemerisProvider> custom_provider_;
};

} // namespace astdyn::propagation

#endif // ASTDYN_HIGH_PRECISION_PROPAGATOR_HPP
