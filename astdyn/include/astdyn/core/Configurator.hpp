#ifndef ASTDYN_CORE_CONFIGURATOR_HPP
#define ASTDYN_CORE_CONFIGURATOR_HPP

#include <string>
#include <istream>
#include <nlohmann/json.hpp>
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/orbit_determination/DifferentialCorrector.hpp"
#include "astdyn/orbit_determination/GaussIOD.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <memory>

namespace astdyn::core {

/**
 * @brief Configuration manager for AstDyn
 * 
 * Capable of parsing standard JSON input to automatically set up 
 * library parameters such as integrators, propagators, and OD configurations.
 */
class Configurator {
public:
    /**
     * @brief Parse JSON to specific settings structures
     */
    static propagation::PropagatorSettings parsePropagatorSettings(const nlohmann::json& j);
    static orbit_determination::DifferentialCorrectorSettings parseDifferentialCorrectorSettings(const nlohmann::json& j);
    static orbit_determination::GaussIODSettings parseGaussIODSettings(const nlohmann::json& j);

    /**
     * @brief Create an Integrator instance from JSON configuration
     */
    static std::shared_ptr<propagation::Integrator> createIntegrator(const nlohmann::json& j);

    /**
     * @brief Create a Propagator instance from JSON configuration
     */
    static std::shared_ptr<propagation::Propagator> createPropagator(const nlohmann::json& j, std::shared_ptr<ephemeris::PlanetaryEphemeris> ephem);
    
    /**
     * @brief Load from standard stream
     */
    static void loadFromStream(std::istream& is,
        propagation::PropagatorSettings& prop_settings,
        orbit_determination::DifferentialCorrectorSettings& dc_settings,
        orbit_determination::GaussIODSettings& iod_settings);

    /**
     * @brief Load from string
     */
    static void loadFromString(const std::string& json_string,
        propagation::PropagatorSettings& prop_settings,
        orbit_determination::DifferentialCorrectorSettings& dc_settings,
        orbit_determination::GaussIODSettings& iod_settings);
};

} // namespace astdyn::core

#endif // ASTDYN_CORE_CONFIGURATOR_HPP
