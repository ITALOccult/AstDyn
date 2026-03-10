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

physics::CartesianStateTyped<core::GCRF> PlanetaryEphemeris::getState(CelestialBody body, time::EpochTDB t) {
    if (body == CelestialBody::SUN) return physics::CartesianStateTyped<core::GCRF>::from_si(
        t,
        math::Vector3<core::GCRF, physics::Distance>::from_si(0,0,0),
        math::Vector3<core::GCRF, physics::Velocity>::from_si(0,0,0)
    );

    if (global_provider_ && global_provider_->isAvailable()) {
        auto sv = global_provider_->getState(body, t);
        auto ss = global_provider_->getState(CelestialBody::SUN, t);
        return physics::CartesianStateTyped<core::GCRF>::from_si(
            t,
            sv.position - ss.position,
            sv.velocity - ss.velocity
        );
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
    return math::Vector3<core::GCRF, physics::Distance>::from_si(0.0, 0.0, 0.0);
}

physics::CartesianStateTyped<core::GCRF> PlanetaryEphemeris::heliocentricToBarycentric(
    const physics::CartesianStateTyped<core::GCRF>& heliocentric_state, 
    time::EpochTDB t) 
{
    auto r_sun = getSunBarycentricPosition(t);
    // getSunBarycentricPosition returns barycentric position of Sun (relative to SSB)
    // bary_pos = helio_pos + r_sun
    return physics::CartesianStateTyped<core::GCRF>::from_si(
        t,
        heliocentric_state.position + r_sun,
        heliocentric_state.velocity // Sun velocity not easily available here, simplified
    );
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
