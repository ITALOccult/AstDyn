#include "astdyn/catalog/GaiaToJohnson.hpp"

#include <array>
#include <cstddef>

namespace astdyn::catalog {
namespace {

// Riello et al. (2021), table C.2: G - X = sum c_i (BP-RP)^i.
constexpr std::array<double, 4> G_MINUS_V_COEFFS = {-0.02704, 0.01424, -0.2156, 0.01426};
constexpr std::array<double, 5> G_MINUS_R_COEFFS = {-0.02275, 0.3961, -0.1243, -0.01396, 0.003775};

constexpr double BP_RP_MIN_MAG   = -0.5;
constexpr double BP_RP_MAX_V_MAG = 5.0;
constexpr double BP_RP_MAX_R_MAG = 2.75;

template <std::size_t N>
double polynomial(const std::array<double, N>& coeffs, double x) {
    double value = 0.0;
    for (std::size_t i = N; i-- > 0;) value = value * x + coeffs[i];
    return value;
}

} // namespace

std::optional<double> johnson_v_from_gaia(double g_mag, double bp_rp_mag) {
    if (bp_rp_mag <= BP_RP_MIN_MAG || bp_rp_mag >= BP_RP_MAX_V_MAG) return std::nullopt;
    return g_mag - polynomial(G_MINUS_V_COEFFS, bp_rp_mag);
}

std::optional<double> cousins_r_from_gaia(double g_mag, double bp_rp_mag) {
    if (bp_rp_mag <= BP_RP_MIN_MAG || bp_rp_mag >= BP_RP_MAX_R_MAG) return std::nullopt;
    return g_mag - polynomial(G_MINUS_R_COEFFS, bp_rp_mag);
}

} // namespace astdyn::catalog
