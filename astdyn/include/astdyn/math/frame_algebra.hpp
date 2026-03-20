#ifndef ASTDYN_MATH_FRAME_ALGEBRA_HPP
#define ASTDYN_MATH_FRAME_ALGEBRA_HPP

#include <Eigen/Core>
#include <cmath>

namespace astdyn::math {

// ============================================================================
// Strongly-Typed Frame Vector
// ============================================================================
template <typename Frame, typename PhysUnit>
class Vector3 {
public:
    // Constructors
    constexpr Vector3() noexcept : x_(0.0), y_(0.0), z_(0.0) {}
    
    // Factory method (expecting standard SI parameters of PhysUnit)
    [[nodiscard]] static constexpr Vector3 from_si(double x, double y, double z) noexcept { 
        return Vector3(x, y, z); 
    }

    // Extraction methods
    [[nodiscard]] constexpr double x_si() const noexcept { return x_; }
    [[nodiscard]] constexpr double y_si() const noexcept { return y_; }
    [[nodiscard]] constexpr double z_si() const noexcept { return z_; }
    
    // Interop with existing Eigen APIs (must acknowledge it returns raw SI values)
    [[nodiscard]] Eigen::Vector3d to_eigen_si() const noexcept { 
        return {x_, y_, z_}; 
    }

    // Magnitudes
    [[nodiscard]] constexpr double squared_norm_si() const noexcept { 
        return x_ * x_ + y_ * y_ + z_ * z_; 
    }
    [[nodiscard]] PhysUnit norm() const noexcept {
        return PhysUnit::from_si(std::sqrt(squared_norm_si()));
    }

    // Arithmetic: Adding two vectors requires SAME Frame and SAME Unit
    [[nodiscard]] constexpr Vector3 operator+(const Vector3& o) const noexcept { 
        return Vector3(x_ + o.x_, y_ + o.y_, z_ + o.z_); 
    }
    [[nodiscard]] constexpr Vector3 operator-(const Vector3& o) const noexcept { 
        return Vector3(x_ - o.x_, y_ - o.y_, z_ - o.z_); 
    }
    constexpr Vector3& operator+=(const Vector3& o) noexcept { 
        x_ += o.x_; y_ += o.y_; z_ += o.z_; return *this; 
    }
    constexpr Vector3& operator-=(const Vector3& o) noexcept { 
        x_ -= o.x_; y_ -= o.y_; z_ -= o.z_; return *this; 
    }
    
    // Scalar operations
    [[nodiscard]] constexpr Vector3 operator*(double scalar) const noexcept { 
        return Vector3(x_ * scalar, y_ * scalar, z_ * scalar); 
    }
    [[nodiscard]] constexpr Vector3 operator/(double scalar) const noexcept { 
        return Vector3(x_ / scalar, y_ / scalar, z_ / scalar); 
    }

    // ========================================================================
    // DELETED operators for mixed Frames to prevent compile-time bugs
    // ========================================================================
    template <typename OtherFrame>
    Vector3<OtherFrame, PhysUnit> operator+(const Vector3<OtherFrame, PhysUnit>&) const = delete;
    
    template <typename OtherFrame>
    Vector3<OtherFrame, PhysUnit> operator-(const Vector3<OtherFrame, PhysUnit>&) const = delete;

private:
    explicit constexpr Vector3(double x, double y, double z) noexcept : x_(x), y_(y), z_(z) {}
    double x_, y_, z_; // Internal representation is always SI units matching PhysUnit
};

// ============================================================================
// Strongly-Typed Rotation Matrix
// ============================================================================
template <typename FromFrame, typename ToFrame>
class RotationMatrix {
public:
    // Factory
    [[nodiscard]] static RotationMatrix from_eigen(Eigen::Matrix3d mat) noexcept { 
        return RotationMatrix(mat); 
    }

    // Rotation Multiplication:
    // Rotates a FromFrame vector, yielding a ToFrame vector with the SAME Unit.
    // The compiler ENFORCES that v is in FromFrame.
    template <typename PhysUnit>
    [[nodiscard]] Vector3<ToFrame, PhysUnit> operator*(const Vector3<FromFrame, PhysUnit>& v) const noexcept {
        const Eigen::Vector3d res = mat_ * v.to_eigen_si();
        return Vector3<ToFrame, PhysUnit>::from_si(res.x(), res.y(), res.z());
    }

    // Adjoint / Inverse: swaps FromFrame and ToFrame
    [[nodiscard]] RotationMatrix<ToFrame, FromFrame> inverse() const noexcept {
        return RotationMatrix<ToFrame, FromFrame>::from_eigen(mat_.transpose());
    }

    // Composition:
    // RotationMatrix<B,C> * RotationMatrix<A,B> -> RotationMatrix<A,C>
    template <typename IntermediateFrame>
    [[nodiscard]] RotationMatrix<IntermediateFrame, ToFrame> 
    operator*(const RotationMatrix<IntermediateFrame, FromFrame>& other) const noexcept {
        return RotationMatrix<IntermediateFrame, ToFrame>::from_eigen(mat_ * other.to_eigen());
    }

    // Extraction
    [[nodiscard]] const Eigen::Matrix3d& to_eigen() const noexcept { return mat_; }

private:
    explicit RotationMatrix(Eigen::Matrix3d mat) noexcept : mat_(mat) {}
    Eigen::Matrix3d mat_;
};

} // namespace astdyn::math

#endif // ASTDYN_MATH_FRAME_ALGEBRA_HPP
