/**
 * @file DE441Provider.hpp
 * @brief JPL DE441 ephemeris provider via CSPICE
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * Requires:
 * - CSPICE library
 * - de441.bsp file (~3.3 GB)
 * 
 * Accuracy: ~cm level
 * Epoch: 1550-2650
 */

#ifndef ASTDYN_DE441_PROVIDER_HPP
#define ASTDYN_DE441_PROVIDER_HPP

#include "EphemerisProvider.hpp"
#include <string>
#include <stdexcept>
#include <memory> 

namespace astdyn::io {
    class SPKReader;
}

namespace astdyn::ephemeris {

/**
 * @brief JPL DE441 ephemeris provider (Native Implementation)
 * 
 * Ultra-precise planetary ephemerides from JPL.
 * Uses native specific SPK reader (Dependency-Free).
 * 
 * Accuracy: ~cm level
 * Speed: Fast (direct access)
 */
class DE441Provider : public EphemerisProvider {
public:
    /**
     * @brief Construct DE441 provider
     * 
     * @param bsp_file Path to de441.bsp file
     * @throws std::runtime_error if file cannot be loaded
     */
    explicit DE441Provider(const std::string& bsp_file);
    
    ~DE441Provider() override;
    
    Eigen::Vector3d getPosition(CelestialBody body, double jd_tdb) override;
    Eigen::Vector3d getVelocity(CelestialBody body, double jd_tdb) override;
    
    std::string getName() const override { return "JPL DE441 (Native)"; }
    double getAccuracy() const override { return 0.001; } // ~1 mas
    bool isAvailable() const override { return loaded_; }
    
private:
    bool loaded_ = false;
    std::string bsp_file_;
    std::unique_ptr<io::SPKReader> reader_;

    /**
     * @brief Convert CelestialBody to NAIF ID
     */
    int bodyToNAIFId(CelestialBody body) const;
    
    /**
     * @brief Convert JD (TDB) to ET (ephemeris time)
     */
    double jdToET(double jd_tdb) const;
    
    /**
     * @brief Rotate from J2000 equatorial to ecliptic
     */
    Eigen::Vector3d equatorialToEcliptic(const Eigen::Vector3d& vec) const;
};

} // namespace astdyn::ephemeris

#endif // ASTDYN_DE441_PROVIDER_HPP
