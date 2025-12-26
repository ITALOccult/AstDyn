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

// Time
#include "astdyn/time/TimeScale.hpp"

// Math
#include "astdyn/math/MathUtils.hpp"

// Ephemeris
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"

// Orbital elements
#include "astdyn/propagation/OrbitalElements.hpp"

// Propagation
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/HighPrecisionPropagator.hpp"

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
