#ifndef ASTDYN_OBSERVATION_MANAGER_HPP
#define ASTDYN_OBSERVATION_MANAGER_HPP

#include "astdyn/observations/Observation.hpp"
#include "astdyn/astrometry/sky_types.hpp"
#include <vector>
#include <string>
#include <map>

namespace astdyn::observations {

/**
 * @brief Manages a collection of astronomical observations.
 * 
 * Handles loading from MPC/RWO formats, sorting by time,
 * and applying catalog bias corrections.
 */
class ObservationManager {
public:
    ObservationManager() = default;

    /**
     * @brief Load observations from MPC format file
     */
    int load_mpc(const std::string& filename);

    /**
     * @brief Add a single observation
     */
    void add(const OpticalObservation& obs);

    /**
     * @brief Clear all observations
     */
    void clear() { observations_.clear(); }

    /**
     * @brief Get reference to observations
     */
    const std::vector<OpticalObservation>& observations() const { return observations_; }

    /**
     * @brief Sort observations by time
     */
    void sort_by_time();

    /**
     * @brief Load catalog biases from CSV (Code, RA_bias, Dec_bias)
     */
    void load_biases(const std::string& filename);

    /**
     * @brief Apply stored biases to all observations
     */
    void apply_biases();

    /**
     * @brief Filter observations by observatory or time range
     */
    void filter_outliers(double sigma_threshold);

private:
    std::vector<OpticalObservation> observations_;
    std::map<std::string, std::pair<astrometry::Angle, astrometry::Angle>> catalog_biases_; // Code -> {RA, Dec} biases
};

} // namespace astdyn::observations

#endif // ASTDYN_OBSERVATION_MANAGER_HPP
