/**
 * @file Vec3.hpp
 * @brief Frame-tagged geometric vectors.
 *
 * Rationale: AstDyn tags *states* by frame (CartesianStateTyped<Frame>) but had
 * no tagged plain vector. Every time the geometry needed a direction or an
 * intermediate position it fell back to a bare Eigen::Vector3d and the frame
 * information evaporated -- which is precisely where the GCRF/CIRS mix-up that
 * cost 16.6 km of ground error was able to compile.
 *
 * Two distinct types, because they are not the same thing:
 *   - Direction<F> : dimensionless unit vector (a line of sight)
 *   - Vec3<F>      : a position or displacement, in metres
 */
#ifndef ASTDYN_MATH_VEC3_HPP
#define ASTDYN_MATH_VEC3_HPP

#include "astdyn/core/frame_tags.hpp"
#include <Eigen/Dense>

namespace astdyn::math {

template <class Frame> class Vec3;

/// Dimensionless unit vector expressed in @c Frame.
template <class Frame>
class Direction {
public:
    Direction() = default;

    /// Normalises on construction; the invariant |d| == 1 always holds.
    [[nodiscard]] static Direction from_eigen(const Eigen::Vector3d& v) {
        return Direction(v.normalized());
    }
    [[nodiscard]] static Direction from_xyz(double x, double y, double z) {
        return from_eigen(Eigen::Vector3d(x, y, z));
    }

    [[nodiscard]] const Eigen::Vector3d& eigen() const noexcept { return v_; }
    [[nodiscard]] double x() const noexcept { return v_.x(); }
    [[nodiscard]] double y() const noexcept { return v_.y(); }
    [[nodiscard]] double z() const noexcept { return v_.z(); }

    /// Dot product is only defined between directions in the SAME frame.
    [[nodiscard]] double dot(const Direction& o) const noexcept { return v_.dot(o.v_); }
    [[nodiscard]] Direction cross(const Direction& o) const { return Direction(v_.cross(o.v_).normalized()); }

    /// Scaling a direction by a length yields a displacement.
    [[nodiscard]] Vec3<Frame> operator*(double metres) const;

private:
    explicit Direction(const Eigen::Vector3d& v) : v_(v) {}
    template <class, class> friend class Rotation;
    Eigen::Vector3d v_{0, 0, 0};
};

/// Position or displacement in @c Frame, stored in metres.
template <class Frame>
class Vec3 {
public:
    Vec3() = default;

    [[nodiscard]] static Vec3 from_m(const Eigen::Vector3d& v) noexcept { return Vec3(v); }
    [[nodiscard]] static Vec3 from_m(double x, double y, double z) noexcept {
        return Vec3(Eigen::Vector3d(x, y, z));
    }
    [[nodiscard]] static Vec3 zero() noexcept { return Vec3(Eigen::Vector3d::Zero()); }

    [[nodiscard]] const Eigen::Vector3d& to_m() const noexcept { return v_; }
    [[nodiscard]] double x_m() const noexcept { return v_.x(); }
    [[nodiscard]] double y_m() const noexcept { return v_.y(); }
    [[nodiscard]] double z_m() const noexcept { return v_.z(); }
    [[nodiscard]] double norm_m() const noexcept { return v_.norm(); }

    /// Same-frame arithmetic only: Vec3<GCRF> + Vec3<ITRF> does not compile.
    [[nodiscard]] Vec3 operator+(const Vec3& o) const noexcept { return Vec3(v_ + o.v_); }
    [[nodiscard]] Vec3 operator-(const Vec3& o) const noexcept { return Vec3(v_ - o.v_); }
    [[nodiscard]] Vec3 operator*(double s) const noexcept { return Vec3(v_ * s); }
    [[nodiscard]] double dot(const Vec3& o) const noexcept { return v_.dot(o.v_); }
    [[nodiscard]] double dot(const Direction<Frame>& d) const noexcept { return v_.dot(d.eigen()); }

    [[nodiscard]] Direction<Frame> direction() const { return Direction<Frame>::from_eigen(v_); }

private:
    explicit Vec3(const Eigen::Vector3d& v) noexcept : v_(v) {}
    template <class, class> friend class Rotation;
    Eigen::Vector3d v_{0, 0, 0};
};

template <class Frame>
Vec3<Frame> Direction<Frame>::operator*(double metres) const {
    return Vec3<Frame>::from_m(v_ * metres);
}

} // namespace astdyn::math

#endif
