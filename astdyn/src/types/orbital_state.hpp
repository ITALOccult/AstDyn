#ifndef ASTDYN_TYPES_ORBAL_STATE_HPP
#define ASTDYN_TYPES_ORBAL_STATE_HPP

#include <array>
#include <type_traits>
#include <optional>

namespace astdyn::types {

// --- Representation Tags ---
struct CartesianTag {};
struct KeplerianTag {};
struct EquinoctialTag {};

/**
 * @brief Immutable template for orbital states.
 * @tparam FrameTag Reference frame tag from core/frame_tags.hpp
 * @tparam RepresentationTag One of the tags defined above.
 */
template <typename FrameTag, typename RepresentationTag>
class OrbitalState {
public:
    using StateData = std::array<double, 6>;

    explicit constexpr OrbitalState(const StateData& data) noexcept : values_(data) {}

    [[nodiscard]] static constexpr std::optional<OrbitalState> create(const StateData& data) noexcept {
        return OrbitalState(data);
    }

    // --- Cartesian Access (requires CartesianTag) ---
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double x() const noexcept { return values_[0]; }
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double y() const noexcept { return values_[1]; }
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double z() const noexcept { return values_[2]; }
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double vx() const noexcept { return values_[3]; }
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double vy() const noexcept { return values_[4]; }
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, CartesianTag>>>
    [[nodiscard]] constexpr double vz() const noexcept { return values_[5]; }

    // --- Keplerian Access (requires KeplerianTag) ---
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double a() const noexcept { return values_[0]; }
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double e() const noexcept { return values_[1]; }
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double i() const noexcept { return values_[2]; }
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double raan() const noexcept { return values_[3]; }
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double arg_peri() const noexcept { return values_[4]; }
    template <typename T = RepresentationTag, typename = std::enable_if_t<std::is_same_v<T, KeplerianTag>>>
    [[nodiscard]] constexpr double m_anomaly() const noexcept { return values_[5]; }

    [[nodiscard]] constexpr const StateData& raw_values() const noexcept { return values_; }

private:
    const StateData values_;
};

} // namespace astdyn::types

#endif // ASTDYN_TYPES_ORBAL_STATE_HPP
