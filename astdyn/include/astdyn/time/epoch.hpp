#ifndef ASTDYN_TIME_EPOCH_HPP
#define ASTDYN_TIME_EPOCH_HPP

#include "astdyn/core/Constants.hpp"

namespace astdyn::time {

// ============================================================================
// Time Scales Tags
// ============================================================================
struct UTC_tag {};
struct TDB_tag {};
struct TT_tag {};
struct TAI_tag {};

// ============================================================================
// Time Duration
// ============================================================================
class TimeDuration {
public:
    // Factory methods
    [[nodiscard]] static constexpr TimeDuration from_seconds(double s) noexcept { return TimeDuration(s); }
    [[nodiscard]] static constexpr TimeDuration from_minutes(double m) noexcept { return TimeDuration(m * 60.0); }
    [[nodiscard]] static constexpr TimeDuration from_hours(double h)   noexcept { return TimeDuration(h * 3600.0); }
    [[nodiscard]] static constexpr TimeDuration from_days(double d)    noexcept { return TimeDuration(d * 86400.0); }
    [[nodiscard]] static constexpr TimeDuration zero()                 noexcept { return TimeDuration(0.0); }

    // Extraction methods
    [[nodiscard]] constexpr double to_seconds() const noexcept { return seconds_; }
    [[nodiscard]] constexpr double to_minutes() const noexcept { return seconds_ / 60.0; }
    [[nodiscard]] constexpr double to_hours()   const noexcept { return seconds_ / 3600.0; }
    [[nodiscard]] constexpr double to_days()    const noexcept { return seconds_ / 86400.0; }

    // Arithmetic operations
    [[nodiscard]] constexpr TimeDuration operator+(const TimeDuration& other) const noexcept { return TimeDuration(seconds_ + other.seconds_); }
    [[nodiscard]] constexpr TimeDuration operator-(const TimeDuration& other) const noexcept { return TimeDuration(seconds_ - other.seconds_); }
    constexpr TimeDuration& operator+=(const TimeDuration& other) noexcept { seconds_ += other.seconds_; return *this; }
    constexpr TimeDuration& operator-=(const TimeDuration& other) noexcept { seconds_ -= other.seconds_; return *this; }
    
    // Scalar operations
    [[nodiscard]] constexpr TimeDuration operator*(double scalar) const noexcept { return TimeDuration(seconds_ * scalar); }
    [[nodiscard]] constexpr TimeDuration operator/(double scalar) const noexcept { return TimeDuration(seconds_ / scalar); }

    // Comparisons
    [[nodiscard]] constexpr bool operator==(const TimeDuration& other) const noexcept { return seconds_ == other.seconds_; }
    [[nodiscard]] constexpr bool operator!=(const TimeDuration& other) const noexcept { return seconds_ != other.seconds_; }
    [[nodiscard]] constexpr bool operator<(const TimeDuration& other) const noexcept { return seconds_ < other.seconds_; }
    [[nodiscard]] constexpr bool operator>(const TimeDuration& other) const noexcept { return seconds_ > other.seconds_; }
    [[nodiscard]] constexpr bool operator<=(const TimeDuration& other) const noexcept { return seconds_ <= other.seconds_; }
    [[nodiscard]] constexpr bool operator>=(const TimeDuration& other) const noexcept { return seconds_ >= other.seconds_; }

private:
    explicit constexpr TimeDuration(double s) noexcept : seconds_(s) {}
    double seconds_; // Internal representation is always seconds (SI)
};

// ============================================================================
// Strongly-Typed Epoch (Time Coordinate)
// ============================================================================
template <typename ScaleTag>
class Epoch {
public:
    // Constructors
    constexpr Epoch() noexcept : mjd_(0.0) {}
    
    // Factory methods
    [[nodiscard]] static constexpr Epoch from_mjd(double mjd_val) noexcept { return Epoch(mjd_val); }
    [[nodiscard]] static constexpr Epoch from_jd(double jd_val)   noexcept { return Epoch(jd_val - 2400000.5); }

    // Extraction methods
    [[nodiscard]] constexpr double mjd() const noexcept { return mjd_; }
    [[nodiscard]] constexpr double jd()  const noexcept { return mjd_ + 2400000.5; }

    // Arithmetic: Delta between two epochs of the SAME scale is a TimeDuration
    [[nodiscard]] constexpr TimeDuration operator-(const Epoch& other) const noexcept {
        return TimeDuration::from_days(mjd_ - other.mjd_);
    }

    // Arithmetic: Advancing an epoch by a TimeDuration yields the same scale
    [[nodiscard]] constexpr Epoch operator+(const TimeDuration& dt) const noexcept {
        return Epoch(mjd_ + dt.to_days());
    }
    [[nodiscard]] constexpr Epoch operator-(const TimeDuration& dt) const noexcept {
        return Epoch(mjd_ - dt.to_days());
    }

    constexpr Epoch& operator+=(const TimeDuration& dt) noexcept { mjd_ += dt.to_days(); return *this; }
    constexpr Epoch& operator-=(const TimeDuration& dt) noexcept { mjd_ -= dt.to_days(); return *this; }

    // Comparisons (only within same scale)
    [[nodiscard]] constexpr bool operator==(const Epoch& other) const noexcept { return mjd_ == other.mjd_; }
    [[nodiscard]] constexpr bool operator!=(const Epoch& other) const noexcept { return mjd_ != other.mjd_; }
    [[nodiscard]] constexpr bool operator<(const Epoch& other) const noexcept { return mjd_ < other.mjd_; }
    [[nodiscard]] constexpr bool operator>(const Epoch& other) const noexcept { return mjd_ > other.mjd_; }
    [[nodiscard]] constexpr bool operator<=(const Epoch& other) const noexcept { return mjd_ <= other.mjd_; }
    [[nodiscard]] constexpr bool operator>=(const Epoch& other) const noexcept { return mjd_ >= other.mjd_; }

    // ========================================================================
    // DELETED operators for mixed scales to prevent compile-time bugs
    // ========================================================================
    template <typename OtherScale>
    TimeDuration operator-(const Epoch<OtherScale>&) const = delete;

    template <typename OtherScale>
    bool operator==(const Epoch<OtherScale>&) const = delete;
    
    template <typename OtherScale>
    bool operator<(const Epoch<OtherScale>&) const = delete;

private:
    explicit constexpr Epoch(double mjd_val) noexcept : mjd_(mjd_val) {}
    double mjd_; // Internal representation is always Modified Julian Date 
};

// ============================================================================
// Type Aliases
// ============================================================================
using EpochUTC = Epoch<UTC_tag>;
using EpochTDB = Epoch<TDB_tag>;
using EpochTT  = Epoch<TT_tag>;
using EpochTAI = Epoch<TAI_tag>;

} // namespace astdyn::time

#endif // ASTDYN_TIME_EPOCH_HPP
