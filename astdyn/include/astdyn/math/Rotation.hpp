/**
 * @file Rotation.hpp
 * @brief Frame-to-frame rotation carrying its endpoints in the type.
 *
 * Supersedes the trait-based sketch in core/frame_rotation.hpp, which could not
 * express the cases we actually need: RotationMatrixProvider<From,To>::apply()
 * is static and therefore has no epoch, while GCRF->ITRF is time dependent.
 * Making the rotation a *value* lets it carry its epoch and compose.
 *
 * The payoff is that the following no longer compiles:
 *
 *     Vec3<GCRF> p = ...;
 *     auto q = era_rotation(jd_ut1) * p;    // ERROR: Rotation<CIRS,TIRS> * Vec3<GCRF>
 *
 * because ERA maps CIRS -> TIRS, not GCRF -> ITRF. Omitting the
 * bias-precession-nutation matrix Q becomes a type error rather than a silent
 * 16.6 km ground error.
 */
#ifndef ASTDYN_MATH_ROTATION_HPP
#define ASTDYN_MATH_ROTATION_HPP

#include "astdyn/math/Vec3.hpp"
#include <Eigen/Dense>

namespace astdyn::math {

template <class From, class To>
class Rotation {
public:
    /**
     * @brief Wrap a raw matrix. Reserved for implementing frame transforms.
     *
     * Deliberately verbose: calling this is asserting, without proof, that @p m
     * really maps @c From to @c To. Every use is a hole in the type system, so
     * they belong in one audited place, not in application code.
     */
    [[nodiscard]] static Rotation from_matrix_unchecked(const Eigen::Matrix3d& m) noexcept {
        return Rotation(m);
    }
    [[nodiscard]] static Rotation identity() noexcept {
        return Rotation(Eigen::Matrix3d::Identity());
    }

    [[nodiscard]] const Eigen::Matrix3d& matrix() const noexcept { return m_; }

    /// Apply to a position: Vec3<From> -> Vec3<To>.
    [[nodiscard]] Vec3<To> operator*(const Vec3<From>& v) const noexcept {
        return Vec3<To>::from_m(m_ * v.to_m());
    }
    /// Apply to a direction: Direction<From> -> Direction<To>.
    [[nodiscard]] Direction<To> operator*(const Direction<From>& d) const {
        return Direction<To>::from_eigen(m_ * d.eigen());
    }

    /// Compose: (Mid->To) * (From->Mid) == (From->To). Only the correct order types.
    template <class Src>
    [[nodiscard]] Rotation<Src, To> operator*(const Rotation<Src, From>& inner) const noexcept {
        return Rotation<Src, To>::from_matrix_unchecked(m_ * inner.matrix());
    }

    /// A rotation is orthogonal, so the inverse is the transpose.
    [[nodiscard]] Rotation<To, From> inverse() const noexcept {
        return Rotation<To, From>::from_matrix_unchecked(m_.transpose());
    }

private:
    explicit Rotation(const Eigen::Matrix3d& m) noexcept : m_(m) {}
    Eigen::Matrix3d m_{Eigen::Matrix3d::Identity()};
};

} // namespace astdyn::math

#endif
