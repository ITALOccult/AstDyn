#pragma once
/**
 * @file GaiaToJohnson.hpp
 * @brief Johnson-Cousins V and R from Gaia (E)DR3 G and BP-RP.
 *
 * Occult4 and OccultWatcher expect Johnson magnitudes in the <Star> record and
 * compute the magnitude drop from them; Gaia G is not V (G - V reaches -0.3 for
 * red stars). The polynomials are the Gaia EDR3 photometric relations of
 * Riello et al. (2021), A&A 649, A3, table C.2, valid only inside their colour
 * range: outside it no value is returned rather than an extrapolation.
 */

#include <optional>

namespace astdyn::catalog {

/// Johnson V from Gaia G and BP-RP; empty outside -0.5 < BP-RP < 5.0.
[[nodiscard]] std::optional<double> johnson_v_from_gaia(double g_mag, double bp_rp_mag);

/// Cousins R from Gaia G and BP-RP; empty outside -0.5 < BP-RP < 2.75.
[[nodiscard]] std::optional<double> cousins_r_from_gaia(double g_mag, double bp_rp_mag);

} // namespace astdyn::catalog
