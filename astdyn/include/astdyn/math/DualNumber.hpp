/**
 * @file DualNumber.hpp
 * @brief Automatic differentiation via Dual Numbers (Vectorized for Gradient).
 */

#ifndef ASTDYN_MATH_DUAL_NUMBER_HPP
#define ASTDYN_MATH_DUAL_NUMBER_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include <array>

namespace astdyn::math {

/**
 * @brief Dual number with N independent variables for gradient computation.
 * 
 * Represents f = val + sum(der[i] * epsilon[i]), where epsilon[i]*epsilon[j] = 0.
 */
template <int N>
class DualV {
public:
    double val;
    std::array<double, N> der;

    constexpr DualV(double v = 0.0) : val(v) {
        der.fill(0.0);
    }

    static constexpr DualV variable(int i, double v) {
        DualV res(v);
        res.der[i] = 1.0;
        return res;
    }

    // Arithmetic
    constexpr DualV operator+(const DualV& o) const {
        DualV res;
        res.val = val + o.val;
        for (int i = 0; i < N; ++i) res.der[i] = der[i] + o.der[i];
        return res;
    }

    constexpr DualV operator-(const DualV& o) const {
        DualV res;
        res.val = val - o.val;
        for (int i = 0; i < N; ++i) res.der[i] = der[i] - o.der[i];
        return res;
    }

    constexpr DualV operator*(const DualV& o) const {
        DualV res;
        res.val = val * o.val;
        for (int i = 0; i < N; ++i) res.der[i] = val * o.der[i] + der[i] * o.val;
        return res;
    }

    constexpr DualV operator/(const DualV& o) const {
        DualV res;
        res.val = val / o.val;
        double o2 = o.val * o.val;
        for (int i = 0; i < N; ++i) res.der[i] = (der[i] * o.val - val * o.der[i]) / o2;
        return res;
    }

    // Scalar arithmetic
    constexpr DualV operator+(double v) const { DualV r = *this; r.val += v; return r; }
    constexpr DualV operator-(double v) const { DualV r = *this; r.val -= v; return r; }
    constexpr DualV operator*(double v) const { DualV r = *this; r.val *= v; for(auto& d : r.der) d *= v; return r; }
    constexpr DualV operator/(double v) const { DualV r = *this; r.val /= v; for(auto& d : r.der) d /= v; return r; }

    friend constexpr DualV operator+(double v, const DualV& d) { return d + v; }
    friend constexpr DualV operator-(double v, const DualV& d) { 
        DualV r; r.val = v - d.val; for(int i=0; i<N; ++i) r.der[i] = -d.der[i]; return r;
    }
    friend constexpr DualV operator*(double v, const DualV& d) { return d * v; }
    friend constexpr DualV operator/(double v, const DualV& d) {
        DualV r; r.val = v / d.val;
        double d2 = d.val * d.val;
        for(int i=0; i<N; ++i) r.der[i] = -v * d.der[i] / d2;
        return r;
    }

    // Math functions
    friend DualV sin(const DualV& d) {
        DualV r; r.val = std::sin(d.val);
        double c = std::cos(d.val);
        for(int i=0; i<N; ++i) r.der[i] = c * d.der[i];
        return r;
    }
    friend DualV cos(const DualV& d) {
        DualV r; r.val = std::cos(d.val);
        double s = -std::sin(d.val);
        for(int i=0; i<N; ++i) r.der[i] = s * d.der[i];
        return r;
    }
    friend DualV sqrt(const DualV& d) {
        DualV r; r.val = std::sqrt(d.val);
        double f = 0.5 / r.val;
        for(int i=0; i<N; ++i) r.der[i] = d.der[i] * f;
        return r;
    }
    friend DualV pow(const DualV& d, double p) {
        DualV r; r.val = std::pow(d.val, p);
        double f = p * std::pow(d.val, p - 1.0);
        for(int i=0; i<N; ++i) r.der[i] = d.der[i] * f;
        return r;
    }

    friend DualV abs(const DualV& d) {
        DualV r; r.val = std::abs(d.val);
        double s = (d.val >= 0) ? 1.0 : -1.0;
        for(int i=0; i<N; ++i) r.der[i] = d.der[i] * s;
        return r;
    }

    friend DualV exp(const DualV& d) {
        DualV r; r.val = std::exp(d.val);
        for(int i=0; i<N; ++i) r.der[i] = r.val * d.der[i];
        return r;
    }

    friend DualV log(const DualV& d) {
        DualV r; r.val = std::log(d.val);
        for(int i=0; i<N; ++i) r.der[i] = d.der[i] / d.val;
        return r;
    }

    friend DualV atan2(const DualV& y, const DualV& x) {
        DualV r; r.val = std::atan2(y.val, x.val);
        double d2 = x.val * x.val + y.val * y.val;
        for(int i=0; i<N; ++i) r.der[i] = (x.val * y.der[i] - y.val * x.der[i]) / d2;
        return r;
    }
};

using Dual6 = DualV<6>;

} // namespace astdyn::math

#endif // ASTDYN_MATH_DUAL_NUMBER_HPP
