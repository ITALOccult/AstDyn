/**
 * @file JPLDEProvider.hpp
 * @brief JPL Development Ephemeris (DE405/DE441) provider using CSPICE
 * 
 * Provides high-accuracy planetary positions from JPL binary ephemeris files.
 * Requires CSPICE toolkit: https://naif.jpl.nasa.gov/naif/toolkit.html
 * 
 * Accuracy:
 * - DE405: ~1 km for planets (1600-2200)
 * - DE441: ~cm level for planets (1550-2650)
 * 
 * Usage:
 *   auto provider = std::make_unique<JPLDEProvider>("de441.bsp");
 *   auto pos = provider->getPosition(CelestialBody::EARTH, jd_tdb);
 */

#ifndef ASTDYN_JPL_DE_PROVIDER_HPP
#define ASTDYN_JPL_DE_PROVIDER_HPP

#include <astdyn/ephemeris/EphemerisInterface.hpp>
#include <string>
#include <map>

// Forward declare CSPICE types to avoid header dependency
extern "C" {
    typedef int SpiceInt;
    typedef double SpiceDouble;
}

namespace astdyn {
namespace ephemeris {

/**
 * @brief JPL DE ephemeris provider using CSPICE
 */
class JPLDEProvider : public IEphemerisProvider {
public:
    /**
     * @brief Construct JPL DE provider
     * @param spk_file Path to SPICE SPK kernel file (e.g., "de441.bsp")
     * @param source Ephemeris source (DE405 or DE441)
     * @throws std::runtime_error if file cannot be loaded
     */
    explicit JPLDEProvider(const std::string& spk_file, 
                           EphemerisSource source = EphemerisSource::JPL_DE441);
    
    /**
     * @brief Destructor - unloads SPICE kernels
     */
    ~JPLDEProvider() override;
    
    // Disable copy (SPICE state is global)
    JPLDEProvider(const JPLDEProvider&) = delete;
    JPLDEProvider& operator=(const JPLDEProvider&) = delete;
    
    Eigen::Vector3d getPosition(CelestialBody body, double jd_tdb) override;
    Eigen::Vector3d getVelocity(CelestialBody body, double jd_tdb) override;
    coordinates::CartesianState getState(CelestialBody body, double jd_tdb) override;
    Eigen::Vector3d getSunBarycentricPosition(double jd_tdb) override;
    
    EphemerisSource getSource() const override { return source_; }
    std::pair<double, double> getValidRange() const override;
    
private:
    /**
     * @brief Convert CelestialBody to NAIF/SPICE body ID
     */
    int getNAIFCode(CelestialBody body) const;
    
    /**
     * @brief Get state vector from SPICE
     * @param target NAIF target body code
     * @param observer NAIF observer body code (10 = Sun)
     * @param et Ephemeris time (TDB seconds past J2000)
     * @return [x, y, z, vx, vy, vz] in km and km/s (J2000 ecliptic)
     */
    std::array<double, 6> getStateFromSPICE(int target, int observer, double et);
    
    /**
     * @brief Convert JD TDB to SPICE ephemeris time (ET)
     * ET = (JD_TDB - 2451545.0) * 86400.0 seconds
     */
    double jdToET(double jd_tdb) const {
        return (jd_tdb - 2451545.0) * 86400.0;
    }
    
    std::string spk_file_;
    EphemerisSource source_;
    int kernel_handle_;  ///< SPICE kernel handle
    
    // NAIF body codes
    static const std::map<CelestialBody, int> naif_codes_;
};

} // namespace ephemeris
} // namespace astdyn

#endif // ASTDYN_JPL_DE_PROVIDER_HPP
