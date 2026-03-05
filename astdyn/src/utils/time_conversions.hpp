#ifndef ASTDYN_UTILS_TIME_CONVERSIONS_HPP
#define ASTDYN_UTILS_TIME_CONVERSIONS_HPP

#include "astdyn/core/Constants.hpp"
#include "time_types.hpp"
#include <optional>

namespace astdyn::utils {

namespace constants = astdyn::constants;

/**
 * @brief Converts TT to TAI.
 * @param tt_instant Instant in TT.
 * @return Instant in TAI, or nullopt if scale mismatch.
 */
[[nodiscard]] constexpr std::optional<Instant> tt_to_tai(const Instant tt_instant) noexcept {
    if (tt_instant.scale != TimeScale::TT) {
        return std::nullopt;
    }

    const double offset_days = constants::TT_MINUS_TAI_SECONDS / constants::SECONDS_PER_DAY;
    const double mjd_tai = tt_instant.mjd.value - offset_days;

    return Instant::from_tai(ModifiedJulianDate(mjd_tai));
}

/**
 * @brief Converts TAI to UTC.
 * @param tai_instant Instant in TAI.
 * @return Instant in UTC, or nullopt if scale mismatch.
 */
[[nodiscard]] constexpr std::optional<Instant> tai_to_utc(const Instant tai_instant) noexcept {
    if (tai_instant.scale != TimeScale::TAI) {
        return std::nullopt;
    }

    const double leap_sec_offset = constants::CURRENT_UTC_TAI_DIFF_SECONDS / constants::SECONDS_PER_DAY;
    const double mjd_utc = tai_instant.mjd.value - leap_sec_offset;

    return Instant::from_utc(ModifiedJulianDate(mjd_utc));
}

/**
 * @brief Direct conversion from TT to UTC.
 * Chain of 15 lines max: Uses internal conversions.
 */
[[nodiscard]] constexpr std::optional<Instant> tt_to_utc(const Instant tt_instant) noexcept {
    const auto tai_opt = tt_to_tai(tt_instant);
    
    if (!tai_opt.has_value()) {
        return std::nullopt;
    }

    return tai_to_utc(tai_opt.value());
}

} // namespace astdyn::utils

#endif // ASTDYN_UTILS_TIME_CONVERSIONS_HPP
