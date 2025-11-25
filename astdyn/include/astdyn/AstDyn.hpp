/**
 * @file AstDyn.hpp
 * @brief Main include file for AstDyn library
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Include this file to access all OrbFit functionality.
 */

#ifndef ASTDYN_HPP
#define ASTDYN_HPP

// Version and configuration
#include "astdyn/Version.hpp"
#include "astdyn/Config.hpp"

// Core types and constants
#include "astdyn/core/Types.hpp"
#include "astdyn/core/Constants.hpp"

// Utilities (to be added in Phase 2)
// #include "astdyn/utils/Logger.hpp"
// #include "astdyn/utils/StringUtils.hpp"

// Time (to be added in Phase 2)
// #include "astdyn/time/TimeScale.hpp"
// #include "astdyn/time/TimeConversions.hpp"

// Math (to be added in Phase 2)
// #include "astdyn/math/MathUtils.hpp"
// #include "astdyn/math/LinearAlgebra.hpp"

// Ephemeris (to be added in Phase 3)
// #include "astdyn/ephemeris/JPLEphemeris.hpp"

// Observations (to be added in Phase 4)
// #include "astdyn/observations/Observation.hpp"
// #include "astdyn/observations/ObservationReader.hpp"

// Orbital elements (to be added in Phase 5)
// #include "astdyn/orbit/KeplerianElements.hpp"
// #include "astdyn/orbit/StateVector.hpp"

// Propagation (to be added in Phase 6)
// #include "astdyn/propagation/OrbitPropagator.hpp"
// #include "astdyn/propagation/ForceModel.hpp"

namespace astdyn {

/**
 * @brief Initialize AstDyn library
 * 
 * Call this function once at program startup to initialize
 * the library and set up any global state.
 * 
 * @return true if initialization was successful
 */
inline bool initialize() {
    // Future: Initialize logging, load configuration files, etc.
    return true;
}

/**
 * @brief Shutdown AstDyn library
 * 
 * Call this function at program exit to clean up resources.
 */
inline void shutdown() {
    // Future: Clean up resources, close files, etc.
}

} // namespace astdyn

#endif // ASTDYN_HPP
