/**
 * @file StellarCovariance.hpp
 * @brief Gaia five-parameter astrometry mapped to the plane of sky (Eq. 17).
 *
 * The shadow position carries the uncertainty of the star as well as that of the
 * asteroid, and the two are independent, so the plane-of-sky covariances add:
 *
 *     C_xi_total = C_xi (dynamics) + C_star
 *
 * The star contribution is the Gaia 5x5 covariance pushed through the map that
 * takes the catalogue parameters to the sky position at the event epoch:
 *
 *     C_star = J_star C5 J_star^T                                      (Eq. 17)
 *
 *     J_star = [ 1  0  P_alpha*  dt   0 ]                              (App. C12)
 *              [ 0  1  P_delta    0  dt ]
 *
 * Two things make this worth its own module. First, the proper-motion term dt
 * dominates: Gaia DR3 is at epoch 2016.0, so an event in 2026 sits 10.6 years
 * away and a star with mu ~ 20 mas/yr has moved 0.2" -- comparable to, or larger
 * than, the closest-approach separation itself. Second, the Gaia covariance is
 * NOT diagonal: the ten correlation coefficients are published precisely because
 * ignoring them misstates the error ellipse.
 *
 * Validated against the NumPy reference in star.py.
 */
#ifndef ASTDYN_ASTROMETRY_STELLAR_COVARIANCE_HPP
#define ASTDYN_ASTROMETRY_STELLAR_COVARIANCE_HPP

#include <Eigen/Dense>
#include <optional>

namespace astdyn::astrometry {

/// Milliarcsecond in radians.
inline constexpr double kMasToRad = 4.84813681109536e-9;
/// Gaia DR3 reference epoch, Julian year.
inline constexpr double kGaiaDR3RefEpoch = 2016.0;

/**
 * @brief The ten Gaia correlation coefficients, in the parameter order
 *        (alpha*, delta, parallax, pmra*, pmdec).
 *
 * Every entry is dimensionless and lives in [-1, 1]. Leaving one unset means
 * "uncorrelated", which is a statement about the star, not a default to reach
 * for lightly: Gaia publishes all ten.
 */
struct GaiaCorrelations {
    std::optional<double> ra_dec;             ///< ra_dec_corr
    std::optional<double> ra_parallax;        ///< ra_parallax_corr
    std::optional<double> ra_pmra;            ///< ra_pmra_corr
    std::optional<double> ra_pmdec;           ///< ra_pmdec_corr
    std::optional<double> dec_parallax;       ///< dec_parallax_corr
    std::optional<double> dec_pmra;           ///< dec_pmra_corr
    std::optional<double> dec_pmdec;          ///< dec_pmdec_corr
    std::optional<double> parallax_pmra;      ///< parallax_pmra_corr
    std::optional<double> parallax_pmdec;     ///< parallax_pmdec_corr
    std::optional<double> pmra_pmdec;         ///< pmra_pmdec_corr
};

/**
 * @brief Dimensionless parallax factors at the event, on the (east, north) basis.
 *
 * Derived from the Earth's barycentric position. Zero neglects the parallax
 * term, which is defensible when the proper-motion term dominates -- but note
 * that a nearby star can have a parallax larger than its ten-year proper motion.
 */
struct ParallaxFactors {
    double p_alpha = 0.0;
    double p_delta = 0.0;
};

/**
 * @brief Assemble the Gaia 5x5 covariance from errors and correlations (Eq. 23).
 *
 *     C_ij = rho_ij sigma_i sigma_j
 *
 * @param errors  (sigma_alpha*, sigma_delta, sigma_parallax) in mas and
 *                (sigma_pmra*, sigma_pmdec) in mas/yr, in that order.
 *                Note that Gaia's ra_error is already in the alpha* sense.
 * @param corr    The ten coefficients; unset entries are taken as zero.
 * @return Covariance in mixed units mas^2 / mas^2 yr^-1 / mas^2 yr^-2.
 * @throws std::invalid_argument on a negative sigma or a correlation outside [-1, 1].
 */
[[nodiscard]] Eigen::Matrix<double, 5, 5> build_c5(
    const Eigen::Matrix<double, 5, 1>& errors,
    const GaiaCorrelations& corr = {});

/// Years between the event epoch and the catalogue reference epoch.
[[nodiscard]] inline double dt_from_epoch(double event_epoch_jyear,
                                          double ref_epoch = kGaiaDR3RefEpoch) {
    return event_epoch_jyear - ref_epoch;
}

/// Star map Jacobian J_star (2x5) at the event epoch (App. C12).
[[nodiscard]] Eigen::Matrix<double, 2, 5> jacobian_star(
    double dt_years, const ParallaxFactors& pf = {});

/**
 * @brief Plane-of-sky stellar covariance at the event epoch (Eq. 17).
 *
 * @param c5         Gaia covariance in mixed mas units, from build_c5().
 * @param dt_years   t_event - ref_epoch.
 * @param pf         Parallax factors.
 * @param to_rad     If true (default) return rad^2, so it adds directly to the
 *                   dynamical C_xi; otherwise mas^2.
 * @return C_star (2x2), order (alpha*, delta).
 */
[[nodiscard]] Eigen::Matrix2d stellar_covariance(
    const Eigen::Matrix<double, 5, 5>& c5,
    double dt_years,
    const ParallaxFactors& pf = {},
    bool to_rad = true);

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_STELLAR_COVARIANCE_HPP
