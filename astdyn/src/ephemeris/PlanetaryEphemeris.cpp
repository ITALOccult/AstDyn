/**
 * @file PlanetaryEphemeris.cpp
 * @brief Implementation of planetary ephemeris calculations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * @warning Internal analytical models (VSOP87/Simon 1994) are DISABLED.
 * This class now requires a high-precision EphemerisProvider (e.g., DE441).
 */

#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/math/MathUtils.hpp"
#include <cmath>
#include "astdyn/coordinates/ReferenceFrame.hpp"

namespace astdyn {
namespace ephemeris {

using Eigen::Vector3d;
using coordinates::CartesianState;

// Convert constants context
using namespace astdyn::core;

// VSOP87/Simon 1994 constants removed (deprecated)

// ============================================================================
// Public Interface
// ============================================================================

std::shared_ptr<EphemerisProvider> PlanetaryEphemeris::global_provider_ = nullptr;

void PlanetaryEphemeris::setProvider(std::shared_ptr<EphemerisProvider> provider) {
    global_provider_ = provider;
}

math::Vector3<core::GCRF, physics::Distance> PlanetaryEphemeris::getPosition(CelestialBody body, time::EpochTDB t) {
    if (body == CelestialBody::SUN) return math::Vector3<core::GCRF, physics::Distance>::from_si(0.0, 0.0, 0.0);

    if (global_provider_ && global_provider_->isAvailable()) {
        auto p_bary = global_provider_->getPosition(body, t);
        auto s_bary = global_provider_->getPosition(CelestialBody::SUN, t);
        return math::Vector3<core::GCRF, physics::Distance>::from_si(
            p_bary.x_si() - s_bary.x_si(),
            p_bary.y_si() - s_bary.y_si(),
            p_bary.z_si() - s_bary.z_si()
        );
    }
    
    throw std::runtime_error("PlanetaryEphemeris: Analytical VSOP87 fallback is DEPRECATED and DISABLED. "
                             "Please use PlanetaryEphemeris::setProvider() to initialize a high-precision source (e.g. DE441).");
}

math::Vector3<core::GCRF, physics::Velocity> PlanetaryEphemeris::getVelocity(CelestialBody body, time::EpochTDB t) {
    if (body == CelestialBody::SUN) return math::Vector3<core::GCRF, physics::Velocity>::from_si(0.0, 0.0, 0.0);

    if (global_provider_ && global_provider_->isAvailable()) {
        auto v_bary = global_provider_->getVelocity(body, t);
        auto vs_bary = global_provider_->getVelocity(CelestialBody::SUN, t);
        return math::Vector3<core::GCRF, physics::Velocity>::from_si(
            v_bary.x_si() - vs_bary.x_si(),
            v_bary.y_si() - vs_bary.y_si(),
            v_bary.z_si() - vs_bary.z_si()
        );
    }
    
    throw std::runtime_error("PlanetaryEphemeris: Analytical VSOP87 fallback is DEPRECATED and DISABLED.");
}

CartesianState PlanetaryEphemeris::getState(CelestialBody body, time::EpochTDB t) {
    if (body == CelestialBody::SUN) return CartesianState(Vector3d::Zero(), Vector3d::Zero());

    if (global_provider_ && global_provider_->isAvailable()) {
        auto sv = global_provider_->getState(body, t);
        auto ss = global_provider_->getState(CelestialBody::SUN, t);
        return CartesianState((sv.head<3>() - ss.head<3>()) / 1000.0, 
                             (sv.tail<3>() - ss.tail<3>()) / 1000.0);
    }

    throw std::runtime_error("PlanetaryEphemeris: no provider set. Call setProvider() first.");
}

math::Vector3<core::GCRF, physics::Distance> PlanetaryEphemeris::getSunBarycentricPosition(time::EpochTDB t) {
    if (global_provider_ && global_provider_->isAvailable()) {
        try {
            return global_provider_->getPosition(CelestialBody::SUN, t);
        } catch (...) {
            return math::Vector3<core::GCRF, physics::Distance>::from_si(0.0, 0.0, 0.0);
        }
    }
    
    // Fallback logic requires analytical planetary positions which are disabled.
    // Return Heliocentric Sun (0,0,0) as best-effort instead of throwing.
    return math::Vector3<core::GCRF, physics::Distance>::from_si(0.0, 0.0, 0.0);
}

CartesianState PlanetaryEphemeris::heliocentricToBarycentric(
    const CartesianState& heliocentric_state, 
    time::EpochTDB t) 
{
    auto r_sun = getSunBarycentricPosition(t);
    
    // Convert heliocentric state (km) back to meters for calculation
    Vector3d r_helio_m = heliocentric_state.position() * 1000.0;
    
    // Barycentric position = heliocentric position - Sun position
    Vector3d r_bary_m = r_helio_m - Vector3d(r_sun.x_si(), r_sun.y_si(), r_sun.z_si());
    
    // Return in km
    return CartesianState(r_bary_m / 1000.0, 
                          heliocentric_state.velocity());
}

// ============================================================================
// Private Implementation
// ============================================================================

// Analytical helper methods removed (deprecated)
void PlanetaryEphemeris::computeOrbitalElements(CelestialBody, double, double[6]) {
    throw std::runtime_error("computeOrbitalElements is deprecated and disabled.");
}

Vector3d PlanetaryEphemeris::elementsToPosition(const double[6]) {
    throw std::runtime_error("elementsToPosition is deprecated and disabled.");
}

Vector3d PlanetaryEphemeris::elementsToVelocity(const double[6], double) {
    throw std::runtime_error("elementsToVelocity is deprecated and disabled.");
}

Vector3d PlanetaryEphemeris::computePerturbations(CelestialBody, double) {
    throw std::runtime_error("computePerturbations is deprecated and disabled.");
}

} // namespace ephemeris
} // namespace astdyn
