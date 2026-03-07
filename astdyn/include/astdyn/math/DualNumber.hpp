/**
 * @file DualNumber.hpp
 * @brief Automatic differentiation via Dual Numbers.
 */

#ifndef ASTDYN_MATH_DUAL_NUMBER_HPP
#define ASTDYN_MATH_DUAL_NUMBER_HPP

#include <cmath>
#include <iostream>

namespace astdyn::math {

/**
 * @brief Simple dual number for automatic differentiation.
 * 
 * Represents f = real + epsilon * dual, where epsilon^2 = 0.
 * In practice, dual is the derivative f'.
 */
class Dual {
public:
    double val;  ///< Function value
    double der;  ///< Derivative value

    Dual(double v = 0.0, double d = 0.0) : val(v), der(d) {}

    // Arithmetic operators
    Dual operator+(const Dual& other) const { return {val + other.val, der + other.der}; }
    Dual operator-(const Dual& other) const { return {val - other.val, der - other.der}; }
    Dual operator*(const Dual& other) const { return {val * other.val, val * other.der + der * other.val}; }
    Dual operator/(const Dual& other) const { return {val / other.val, (der * other.val - val * other.der) / (other.val * other.val)}; }

    // Constants
    Dual operator+(double v) const { return {val + v, der}; }
    Dual operator-(double v) const { return {val - v, der}; }
    Dual operator*(double v) const { return {val * v, der * v}; }
    Dual operator/(double v) const { return {val / v, der / v}; }

    friend Dual operator+(double v, const Dual& d) { return d + v; }
    friend Dual operator-(double v, const Dual& d) { return {v - d.val, -d.der}; }
    friend Dual operator*(double v, const Dual& d) { return d * v; }
    friend Dual operator/(double v, const Dual& d) { return {v / d.val, -v * d.der / (d.val * d.val)}; }

    // Math functions
    friend Dual sin(const Dual& d) { return {std::sin(d.val), std::cos(d.val) * d.der}; }
    friend Dual cos(const Dual& d) { return {std::cos(d.val), -std::sin(d.val) * d.der}; }
    friend Dual sqrt(const Dual& d) { double s = std::sqrt(d.val); return {s, 0.5 * d.der / s}; }
    friend Dual atan2(const Dual& y, const Dual& x) {
        double r2 = x.val * x.val + y.val * y.val;
        return {std::atan2(y.val, x.val), (y.der * x.val - y.val * x.der) / r2};
    }
    friend Dual asin(const Dual& d) { return {std::asin(d.val), d.der / std::sqrt(1.0 - d.val * d.val)}; }
    
};

} // namespace astdyn::math

#endif // ASTDYN_MATH_DUAL_NUMBER_HPP
