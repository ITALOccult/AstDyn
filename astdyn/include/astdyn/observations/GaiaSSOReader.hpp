/**
 * @file GaiaSSOReader.hpp
 * @brief Gaia asteroid astrometry (DR3 / FPR table sso_observation, CSV).
 */
#ifndef ASTDYN_OBSERVATIONS_GAIA_SSO_READER_HPP
#define ASTDYN_OBSERVATIONS_GAIA_SSO_READER_HPP

#include "astdyn/observations/Observation.hpp"
#include <string>
#include <vector>

namespace astdyn::observations {

/**
 * @brief Reads the CSV of gaiafpr.sso_observation (or gaiadr3.) as returned by
 *        the Gaia archive TAP service, columns addressed by name.
 *
 * These observations are not distributed by the MPC nor by AstDyS. On (3830)
 * the 626 of the FPR are what JPL has and we did not: JPL's orbit fits them to
 * 4 mas along the scan, AstDyS's to 37 mas.
 *
 * Each record carries Gaia's geocentric position (au) and the full
 * (alpha*, delta) covariance, systematic plus random -- the error ellipse is a
 * segment along the scan, and only the correlated covariance keeps that
 * precision. Records flagged is_rejected by the Gaia pipeline are skipped.
 */
std::vector<OpticalObservation> read_gaia_sso_csv(const std::string& path);

} // namespace astdyn::observations

#endif // ASTDYN_OBSERVATIONS_GAIA_SSO_READER_HPP
