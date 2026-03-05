#ifndef ASTDYN_CORE_UNITS_HPP
#define ASTDYN_CORE_UNITS_HPP

namespace astdyn::core {

/**
 * @brief Strong type for Meter.
 */
struct Meter {
    const double value;
    explicit constexpr Meter(const double v) : value(v) {}
};

/**
 * @brief Strong type for Second.
 */
struct Second {
    const double value;
    explicit constexpr Second(const double v) : value(v) {}
};

/**
 * @brief Strong type for Radian.
 */
struct Radian {
    const double value;
    explicit constexpr Radian(const double v) : value(v) {}
};

/**
 * @brief Strong type for Degree.
 */
struct Degree {
    const double value;
    explicit constexpr Degree(const double v) : value(v) {}
};

/**
 * @brief Strong type for MilliArcSecond.
 */
struct MilliArcSecond {
    const double value;
    explicit constexpr MilliArcSecond(const double v) : value(v) {}
};

} // namespace astdyn::core

#endif // ASTDYN_CORE_UNITS_HPP
