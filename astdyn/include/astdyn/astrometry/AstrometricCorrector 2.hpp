#ifndef ASTDYN_ASTROMETRY_ASTROMETRIC_CORRECTOR_HPP
#define ASTDYN_ASTROMETRY_ASTROMETRIC_CORRECTOR_HPP

#include "astdyn/astrometry/AstrometricCorrections.hpp"
#include "astdyn/astrometry/AstrometricTypes.hpp"

namespace astdyn::astrometry {

/**
 * @brief High-level corrector to apply physical astrometric effects (CTFIYH).
 * 
 * This class follows the principle of "Code that fits in your head" by providing
 * a single entry point for corrections, hiding the complexity of frame transformations
 * and sequence of application.
 */
class AstrometricCorrector {
public:
    explicit AstrometricCorrector(const AstrometricSettings& settings)
        : settings_(settings) {}

    /**
     * @brief Apply all enabled corrections to a geocentric equatorial vector.
     * 
     * @param rho_eq_icrf Geocentric vector in Equatorial J2000/ICRF (meters).
     * @param earth_vel_eq Earth velocity (m/s).
     * @param obs_to_sun_eq Vector from observer to Sun (meters).
     * @return Corrected vector.
     */
    Eigen::Vector3d apply(
        const Eigen::Vector3d& rho_eq_icrf,
        const Eigen::Vector3d& earth_vel_eq,
        const Eigen::Vector3d& obs_to_sun_eq) const 
    {
        Eigen::Vector3d corrected = rho_eq_icrf;
        
        // 1. Relativistic Light Deflection (Apply to geometric vector)
        if (settings_.deflessione_relativistica) {
            corrected = deflessione_relativistica(corrected, obs_to_sun_eq);
        }
        
        // 2. Differential Aberration (Apply to deflected vector)
        if (settings_.aberrazione_differenziale) {
            corrected = aberrazione_differenziale(corrected, earth_vel_eq);
        }
        
        return corrected;
    }

private:
    AstrometricSettings settings_;
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_ASTROMETRIC_CORRECTOR_HPP
