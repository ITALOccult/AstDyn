#ifndef ASTDYN_PROPAGATION_FORCE_FIELD_HPP
#define ASTDYN_PROPAGATION_FORCE_FIELD_HPP

#include "astdyn/core/Constants.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/AsteroidPerturbations.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include <Eigen/Core>
#include <memory>
#include <vector>

namespace astdyn::propagation {

/**
 * @brief Base interface for force model components (CTFYH Standard)
 */
class ForceModel {
public:
    virtual ~ForceModel() = default;
    
    /**
     * @brief Compute acceleration [AU/day^2]
     * @param t Time (TDB)
     * @param pos_au Particle position [AU] (Heliocentric or Barycentric depending on field)
     * @param vel_au_d Particle velocity [AU/day]
     * @return Acceleration [AU/day^2]
     */
    virtual Eigen::Vector3d compute_acceleration(
        time::EpochTDB t, 
        const Eigen::Vector3d& pos_au, 
        const Eigen::Vector3d& vel_au_d) const = 0;
};

/**
 * @brief Standard Force Field aggregator
 */
class ForceField {
public:
    explicit ForceField(const PropagatorSettings& settings, 
                        std::shared_ptr<ephemeris::PlanetaryEphemeris> ephem = nullptr);

    /**
     * @brief Compute total acceleration including all active models
     * Output is in AU/day^2.
     */
    Eigen::Vector3d total_acceleration(
        time::EpochTDB t, 
        const Eigen::Vector3d& pos_au, 
        const Eigen::Vector3d& vel_au_d) const;

    const PropagatorSettings& settings() const { return settings_; }
    void update_settings(const PropagatorSettings& s) { settings_ = s; }

private:
    PropagatorSettings settings_;
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
    std::shared_ptr<ephemeris::AsteroidPerturbations> asteroids_;

    // Individual force calculations (internal refactoring for reuse)
    Eigen::Vector3d n_body_perturbation(time::EpochTDB t, const Eigen::Vector3d& pos_au) const;
    Eigen::Vector3d relativistic_correction(const Eigen::Vector3d& r, const Eigen::Vector3d& v) const;
    Eigen::Vector3d j2_correction(time::EpochTDB t, const Eigen::Vector3d& pos_au) const;
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATION_FORCE_FIELD_HPP
