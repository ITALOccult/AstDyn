#ifndef ASTDYN_ASTROMETRY_ASTROMETRIC_CORRECTIONS_HPP
#define ASTDYN_ASTROMETRY_ASTROMETRIC_CORRECTIONS_HPP

#include <Eigen/Dense>

namespace astdyn::astrometry {

/**
 * @brief Applica l'aberrazione differenziale (Stellar/Planetary Aberration).
 * 
 * @param rho_eq Vettore posizione geocentrico (Equatoriale J2000).
 * @param earth_velocity_eq Vettore velocità della Terra (Equatoriale J2000, km/s).
 * @return Vettore corretto per l'aberrazione.
 */
Eigen::Vector3d aberrazione_differenziale(
    const Eigen::Vector3d& rho_eq,
    const Eigen::Vector3d& earth_velocity_eq);

/**
 * @brief Applica la deflessione relativistica (Gravitational Light Deflection).
 * 
 * @param rho_eq Vettore posizione geocentrico (Equatoriale J2000).
 * @param observer_to_sun_eq Vettore dall'osservatore al Sole (Equatoriale J2000).
 * @return Vettore corretto per la deflessione gravitazionale.
 */
Eigen::Vector3d deflessione_relativistica(
    const Eigen::Vector3d& rho_eq,
    const Eigen::Vector3d& observer_to_sun_eq);

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_ASTROMETRIC_CORRECTIONS_HPP
