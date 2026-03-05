#ifndef ASTDYN_UTILS_TIME_TYPES_HPP
#define ASTDYN_UTILS_TIME_TYPES_HPP

#include <concepts>

namespace astdyn::utils {

/**
 * @brief Strong type for Modified Julian Date (MJD).
 */
struct ModifiedJulianDate {
    double value;
    constexpr ModifiedJulianDate() noexcept : value(0.0) {}
    explicit constexpr ModifiedJulianDate(const double val) noexcept : value(val) {}
};

/**
 * @brief Strong type for Julian Date (JD).
 */
struct JulianDate {
    double value;
    constexpr JulianDate() noexcept : value(0.0) {}
    explicit constexpr JulianDate(const double val) noexcept : value(val) {}
};

/**
 * @brief Enumeration of fundamental time scales.
 */
enum class TimeScale {
    UTC, ///< Universal Time Coordinated
    TAI, ///< Atomic Time
    TT,  ///< Terrestrial Time
    UT1  ///< Universal Time 1
};

/**
 * @brief Representation of a temporal point in a specific scale.
 * 
 * CTFYH adherence: Private constructor, Factory methods for data honesty.
 */
struct Instant {
    ModifiedJulianDate mjd;
    TimeScale scale;

    /** @brief Default constructor: J2000.0 (TT) */
    constexpr Instant() noexcept : mjd(51544.5), scale(TimeScale::TT) {}

    /** @brief Factory for UTC instantiation. */
    [[nodiscard]] static constexpr Instant from_utc(const ModifiedJulianDate mjd) noexcept {
        return Instant(mjd, TimeScale::UTC);
    }

    /** @brief Factory for TAI instantiation. */
    [[nodiscard]] static constexpr Instant from_tai(const ModifiedJulianDate mjd) noexcept {
        return Instant(mjd, TimeScale::TAI);
    }

    /** @brief Factory for TT instantiation. */
    [[nodiscard]] static constexpr Instant from_tt(const ModifiedJulianDate mjd) noexcept {
        return Instant(mjd, TimeScale::TT);
    }

    /** @brief Factory for UT1 instantiation. */
    [[nodiscard]] static constexpr Instant from_ut1(const ModifiedJulianDate mjd) noexcept {
        return Instant(mjd, TimeScale::UT1);
    }

private:
    explicit constexpr Instant(const ModifiedJulianDate mjd_val, const TimeScale scale_val) noexcept
        : mjd(mjd_val), scale(scale_val) {}
};

} // namespace astdyn::utils

#endif // ASTDYN_UTILS_TIME_TYPES_HPP
