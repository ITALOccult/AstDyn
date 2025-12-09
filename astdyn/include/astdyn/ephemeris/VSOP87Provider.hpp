/**
 * @file VSOP87Provider.hpp
 * @brief VSOP87 ephemeris provider (built-in, no dependencies)
 * @author AstDyn Team
 * @date 2025-12-09
 */

#ifndef ASTDYN_VSOP87_PROVIDER_HPP
#define ASTDYN_VSOP87_PROVIDER_HPP

#include "EphemerisProvider.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"

namespace astdyn::ephemeris {

/**
 * @brief VSOP87 analytical ephemeris provider
 * 
 * Built-in provider using VSOP87 analytical series.
 * No external dependencies, fast, good accuracy for 1800-2050.
 * 
 * Accuracy: ~1-20 arcsec
 * Speed: Very fast
 * Dependencies: None
 */
class VSOP87Provider : public EphemerisProvider {
public:
    VSOP87Provider() = default;
    
    Eigen::Vector3d getPosition(CelestialBody body, double jd_tdb) override {
        return PlanetaryEphemeris::getPosition(body, jd_tdb);
    }
    
    Eigen::Vector3d getVelocity(CelestialBody body, double jd_tdb) override {
        return PlanetaryEphemeris::getVelocity(body, jd_tdb);
    }
    
    std::string getName() const override {
        return "VSOP87";
    }
    
    double getAccuracy() const override {
        return 20.0;  // ~20 arcsec
    }
    
    bool isAvailable() const override {
        return true;  // Always available (built-in)
    }
};

} // namespace astdyn::ephemeris

#endif // ASTDYN_VSOP87_PROVIDER_HPP
