#ifndef ASTDYN_CORE_UNITS_HPP
#define ASTDYN_CORE_UNITS_HPP

namespace astdyn::core {

/**
 * @brief Strong type for Meter.
 */
struct Meter {
    double value;
    explicit constexpr Meter(const double v) : value(v) {}
    constexpr Meter() : value(0.0) {}
};

/**
 * @brief Strong type for Second.
 */
struct Second {
    double value;
    explicit constexpr Second(const double v) : value(v) {}
    constexpr Second() : value(0.0) {}
};

/**
 * @brief Strong type for Radian.
 */
struct Radian {
    double value;
    explicit constexpr Radian(const double v) : value(v) {}
    constexpr Radian() : value(0.0) {}
};

/**
 * @brief Strong type for Degree.
 */
struct Degree {
    double value;
    explicit constexpr Degree(const double v) : value(v) {}
    constexpr Degree() : value(0.0) {}
};

/**
 * @brief Strong type for MilliArcSecond.
 */
struct MilliArcSecond {
    double value;
    explicit constexpr MilliArcSecond(const double v) : value(v) {}
    constexpr MilliArcSecond() : value(0.0) {}
};

} // namespace astdyn::core

#endif // ASTDYN_CORE_UNITS_HPP
