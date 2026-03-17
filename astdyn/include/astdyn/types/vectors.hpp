#ifndef ASTDYN_TYPES_VECTORS_HPP
#define ASTDYN_TYPES_VECTORS_HPP

#include <Eigen/Dense>
#include <cmath>

namespace astdyn::types {

/**
 * @brief Strong-typed 3D Vector with Reference Frame and Physical Unit tagging.
 * @tparam Frame Reference frame tag (e.g., core::GCRF, core::ITRF).
 * @tparam Unit Physical unit tag (e.g., core::Meter).
 */
template <typename Frame, typename Unit>
struct Vector3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr Vector3() noexcept = default;
    explicit constexpr Vector3(const double v1, const double v2, const double v3) noexcept 
        : x(v1), y(v2), z(v3) {}

    /** @brief Construct from Eigen vector. */
    explicit Vector3(const Eigen::Vector3d& v) noexcept 
        : x(v.x()), y(v.y()), z(v.z()) {}

    /** @brief Convert to Eigen vector. */
    [[nodiscard]] Eigen::Vector3d to_eigen() const noexcept {
        return Eigen::Vector3d(x, y, z);
    }

    /** @brief Dot product. */
    [[nodiscard]] constexpr double dot(const Vector3& other) const noexcept {
        return (x * other.x) + (y * other.y) + (z * other.z);
    }

    /** @brief Cross product. */
    [[nodiscard]] constexpr Vector3 cross(const Vector3& other) const noexcept {
        return Vector3(
            (y * other.z) - (z * other.y),
            (z * other.x) - (x * other.z),
            (x * other.y) - (y * other.x)
        );
    }

    /** @brief Vector norm (magnitude). */
    [[nodiscard]] double norm() const noexcept {
        return std::sqrt(x * x + y * y + z * z);
    }

    /** @brief Squared norm. */
    [[nodiscard]] constexpr double squaredNorm() const noexcept {
        return x * x + y * y + z * z;
    }

    /** @brief Normalized vector (unit vector). */
    [[nodiscard]] Vector3 normalized() const noexcept {
        double n = norm();
        if (n < 1e-18) return Vector3(0, 0, 0);
        return Vector3(x / n, y / n, z / n);
    }

    /** @brief Subtraction. */
    [[nodiscard]] constexpr Vector3 operator-(const Vector3& other) const noexcept {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    /** @brief Addition. */
    [[nodiscard]] constexpr Vector3 operator+(const Vector3& other) const noexcept {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    /** @brief Negation. */
    [[nodiscard]] constexpr Vector3 operator-() const noexcept {
        return Vector3(-x, -y, -z);
    }

    /** @brief Scalar multiplication. */
    [[nodiscard]] constexpr Vector3 operator*(double scalar) const noexcept {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }

    /** @brief Scalar division. */
    [[nodiscard]] constexpr Vector3 operator/(double scalar) const noexcept {
        return Vector3(x / scalar, y / scalar, z / scalar);
    }

    /** @brief In-place addition. */
    constexpr Vector3& operator+=(const Vector3& other) noexcept {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    /** @brief In-place subtraction. */
    constexpr Vector3& operator-=(const Vector3& other) noexcept {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    /** @brief In-place addition from Eigen vector (component-wise). */
    Vector3& operator+=(const Eigen::Vector3d& other) noexcept {
        x += other.x();
        y += other.y();
        z += other.z();
        return *this;
    }
    
    /** @brief Friend scalar multiplication (scalar * vector). */
    friend constexpr Vector3 operator*(double scalar, const Vector3& vec) noexcept {
        return vec * scalar;
    }
};

} // namespace astdyn::types

#endif // ASTDYN_TYPES_VECTORS_HPP
