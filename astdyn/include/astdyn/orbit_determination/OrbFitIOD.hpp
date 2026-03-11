#ifndef ASTDYN_ORBIT_DETERMINATION_ORBFIT_IOD_HPP
#define ASTDYN_ORBIT_DETERMINATION_ORBFIT_IOD_HPP

#include "astdyn/orbit_determination/GaussIOD.hpp"
#include <vector>

namespace astdyn::orbit_determination {

/**
 * @brief OrbFit's robust Gauss/Laplace Initial Orbit Determination methods.
 * 
 * Provides AstDyn-compatible wrappers around the C++ OrbFit analytical solutions.
 */
class OrbFitIOD {
public:
    explicit OrbFitIOD(const GaussIODSettings& settings = GaussIODSettings());

    /**
     * @brief Determine orbit from multiple observations (selects best 3).
     * Uses OrbFit Gauss method by default, falling back to Laplace if needed.
     */
    GaussIODResult compute(const std::vector<observations::OpticalObservation>& observations);

    /**
     * @brief Determine orbit exactly from three observations.
     * 
     * @param obs1 First observation
     * @param obs2 Second observation (epoch of solution)
     * @param obs3 Third observation
     */
    GaussIODResult compute_from_three(
        const observations::OpticalObservation& obs1,
        const observations::OpticalObservation& obs2,
        const observations::OpticalObservation& obs3);

    const GaussIODSettings& settings() const { return settings_; }
    void set_settings(const GaussIODSettings& s) { settings_ = s; }

private:
    GaussIODSettings settings_;

    std::optional<std::array<int, 3>> select_observations(
        const std::vector<observations::OpticalObservation>& observations) const;
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_ORBFIT_IOD_HPP
