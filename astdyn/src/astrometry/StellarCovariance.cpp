#include "astdyn/astrometry/StellarCovariance.hpp"
#include <stdexcept>
#include <string>
#include <array>
#include <utility>

namespace astdyn::astrometry {
namespace {

/// Check and apply one correlation coefficient.
void set_corr(Eigen::Matrix<double, 5, 5>& C, const Eigen::Matrix<double, 5, 1>& sig,
              int i, int j, const std::optional<double>& rho, const char* name) {
    if (!rho) return;
    if (*rho < -1.0 || *rho > 1.0) {
        throw std::invalid_argument(std::string("StellarCovariance: ") + name + " = "
                                    + std::to_string(*rho) + " outside [-1, 1]");
    }
    C(i, j) = C(j, i) = *rho * sig(i) * sig(j);
}

} // namespace

Eigen::Matrix<double, 5, 5> build_c5(const Eigen::Matrix<double, 5, 1>& errors,
                                     const GaiaCorrelations& corr) {
    if ((errors.array() < 0.0).any()) {
        throw std::invalid_argument("StellarCovariance: sigmas must be non-negative");
    }
    Eigen::Matrix<double, 5, 5> C = Eigen::Matrix<double, 5, 5>::Zero();
    C.diagonal() = errors.array().square();

    set_corr(C, errors, 0, 1, corr.ra_dec,          "ra_dec_corr");
    set_corr(C, errors, 0, 2, corr.ra_parallax,     "ra_parallax_corr");
    set_corr(C, errors, 0, 3, corr.ra_pmra,         "ra_pmra_corr");
    set_corr(C, errors, 0, 4, corr.ra_pmdec,        "ra_pmdec_corr");
    set_corr(C, errors, 1, 2, corr.dec_parallax,    "dec_parallax_corr");
    set_corr(C, errors, 1, 3, corr.dec_pmra,        "dec_pmra_corr");
    set_corr(C, errors, 1, 4, corr.dec_pmdec,       "dec_pmdec_corr");
    set_corr(C, errors, 2, 3, corr.parallax_pmra,   "parallax_pmra_corr");
    set_corr(C, errors, 2, 4, corr.parallax_pmdec,  "parallax_pmdec_corr");
    set_corr(C, errors, 3, 4, corr.pmra_pmdec,      "pmra_pmdec_corr");
    return C;
}

Eigen::Matrix<double, 2, 5> jacobian_star(double dt_years, const ParallaxFactors& pf) {
    Eigen::Matrix<double, 2, 5> J = Eigen::Matrix<double, 2, 5>::Zero();
    J(0, 0) = 1.0;
    J(0, 2) = pf.p_alpha;
    J(0, 3) = dt_years;
    J(1, 1) = 1.0;
    J(1, 2) = pf.p_delta;
    J(1, 4) = dt_years;
    return J;
}

Eigen::Matrix2d stellar_covariance(const Eigen::Matrix<double, 5, 5>& c5,
                                   double dt_years, const ParallaxFactors& pf,
                                   bool to_rad) {
    const auto J = jacobian_star(dt_years, pf);
    Eigen::Matrix2d C = J * c5 * J.transpose();          // mas^2
    if (to_rad) C *= kMasToRad * kMasToRad;              // rad^2
    return 0.5 * (C + C.transpose());                    // symmetrise against round-off
}

} // namespace astdyn::astrometry
