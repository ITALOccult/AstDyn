#ifndef ASTDYN_TYPES_ORBAL_STATE_HPP
#define ASTDYN_TYPES_ORBAL_STATE_HPP

#include <array>
#include <type_traits>
#include <optional>
#include "astdyn/core/Constants.hpp"

namespace astdyn::types {

// --- Representation Tags ---
struct CartesianTag {};
struct KeplerianTag {};
struct EquinoctialTag {};

// --- Unit Tags (for Cartesian states) ---
// Position in meters, velocity in m/s
struct MeterTag {};
// Position in km, velocity in km/s  
struct KilometerTag {};
// Position in AU, velocity in AU/day
struct AuDayTag {};

/**
 * @brief Immutable strongly-typed orbital state.
 *
 * @tparam FrameTag       Reference frame tag (e.g. core::GCRF, core::ECLIPJ2000).
 * @tparam RepTag         Representation tag: CartesianTag, KeplerianTag, EquinoctialTag.
 * @tparam UnitTag        Unit tag for Cartesian states: MeterTag (default), KilometerTag, AuDayTag.
 *                        Ignored for Keplerian (elements always in AU/rad).
 *
 * The UnitTag parameter makes the physical units part of the type, so that
 * passing km-valued vectors to a function expecting meters is a compile-time error.
 */
template <typename FrameTag, typename RepTag, typename UnitTag = MeterTag>
class OrbitalState {
public:
    using StateData = std::array<double, 6>;

    explicit constexpr OrbitalState(const StateData& data) noexcept : values_(data) {}

    [[nodiscard]] static constexpr std::optional<OrbitalState> create(const StateData& data) noexcept {
        return OrbitalState(data);
    }

    // --- Cartesian Access (requires CartesianTag) ---
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double x() const noexcept { return values_[0]; }
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double y() const noexcept { return values_[1]; }
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double z() const noexcept { return values_[2]; }
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double vx() const noexcept { return values_[3]; }
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double vy() const noexcept { return values_[4]; }
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double vz() const noexcept { return values_[5]; }

    // --- Keplerian Access (requires KeplerianTag) ---
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double a() const noexcept { return values_[0]; }
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double e() const noexcept { return values_[1]; }
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double i() const noexcept { return values_[2]; }
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double raan() const noexcept { return values_[3]; }
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double arg_peri() const noexcept { return values_[4]; }
    template <typename T = RepTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double m_anomaly() const noexcept { return values_[5]; }

    [[nodiscard]] constexpr const StateData& raw_values() const noexcept { return values_; }

    /**
     * @brief Converts this Cartesian state to meters/m·s.
     * Returns a new OrbitalState<FrameTag, CartesianTag, MeterTag>.
     * Only meaningful for CartesianTag states.
     */
    template <typename T = RepTag, typename U = UnitTag,
              typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] OrbitalState<FrameTag, CartesianTag, MeterTag> to_meters() const noexcept {
        if constexpr (std::is_same_v<U, MeterTag>) {
            return OrbitalState<FrameTag, CartesianTag, MeterTag>(values_);
        } else if constexpr (std::is_same_v<U, KilometerTag>) {
            return OrbitalState<FrameTag, CartesianTag, MeterTag>({
                values_[0] * 1000.0, values_[1] * 1000.0, values_[2] * 1000.0,
                values_[3] * 1000.0, values_[4] * 1000.0, values_[5] * 1000.0
            });
        } else { // AuDayTag
            constexpr double AU_M = constants::AU * 1000.0;
            constexpr double AU_M_S = AU_M / constants::SECONDS_PER_DAY;
            return OrbitalState<FrameTag, CartesianTag, MeterTag>({
                values_[0] * AU_M,   values_[1] * AU_M,   values_[2] * AU_M,
                values_[3] * AU_M_S, values_[4] * AU_M_S, values_[5] * AU_M_S
            });
        }
    }

private:
    StateData values_;
};

} // namespace astdyn::types

#endif // ASTDYN_TYPES_ORBAL_STATE_HPP
