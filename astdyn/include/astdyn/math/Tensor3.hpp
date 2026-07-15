/**
 * @file Tensor3.hpp
 * @brief Small dense rank-3 tensors for second-order variational calculus.
 *
 * Tensor3 holds the spatial third derivative of the potential, U_ijk
 * (3x3x3, fully symmetric).  Tensor6 holds the state-transition tensor
 * (STT) Psi_iab (6x6x6).  Both are thin wrappers over Eigen fixed-size
 * matrices, following the AstDyn convention of Eigen-in-internals.
 */

#ifndef ASTDYN_MATH_TENSOR3_HPP
#define ASTDYN_MATH_TENSOR3_HPP

#include "astdyn/core/Types.hpp"
#include <Eigen/Dense>
#include <array>

namespace astdyn::math {

/**
 * @brief Fully symmetric rank-3 spatial tensor T_ijk (i,j,k in {0,1,2}).
 *
 * Stored slice-wise: slice(i)(j,k) = T_ijk.  The bilinear contraction
 * T_ijk B_ja B_kb (used by the STT Kick source) is a per-slice
 * congruence B^T slice(i) B and is provided directly.
 */
class Tensor3 {
public:
    Tensor3() { setZero(); }

    void setZero() { for (auto& s : slices_) s.setZero(); }

    double  operator()(int i, int j, int k) const { return slices_[i](j, k); }
    double& operator()(int i, int j, int k)       { return slices_[i](j, k); }

    const astdyn::Matrix3d& slice(int i) const { return slices_[i]; }
    astdyn::Matrix3d&       slice(int i)       { return slices_[i]; }

    Tensor3& operator+=(const Tensor3& o) {
        for (int i = 0; i < 3; ++i) slices_[i] += o.slices_[i];
        return *this;
    }

    Tensor3& operator*=(double s) {
        for (auto& m : slices_) m *= s;
        return *this;
    }

    /**
     * @brief Congruence of slice i with a 3xN block: (B^T slice(i) B).
     * With B = position rows of the STM (3x6), this returns the 6x6
     * source matrix S_i(a,b) = T_ijk B_ja B_kb for momentum row i.
     */
    template <int N>
    Eigen::Matrix<double, N, N>
    contract_bilinear(int i, const Eigen::Matrix<double, 3, N>& B) const {
        return B.transpose() * slices_[i] * B;
    }

private:
    std::array<astdyn::Matrix3d, 3> slices_;   ///< slice(i)(j,k) = T_ijk
};

/**
 * @brief State-transition tensor Psi_iab (6x6x6), stored slice-wise:
 * slice(i)(a,b) = Psi_iab.
 */
class Tensor6 {
public:
    Tensor6() { setZero(); }

    void setZero() { for (auto& s : slices_) s.setZero(); }

    double  operator()(int i, int a, int b) const { return slices_[i](a, b); }
    double& operator()(int i, int a, int b)       { return slices_[i](a, b); }

    const astdyn::Matrix6d& slice(int i) const { return slices_[i]; }
    astdyn::Matrix6d&       slice(int i)       { return slices_[i]; }

private:
    std::array<astdyn::Matrix6d, 6> slices_;   ///< slice(i)(a,b) = Psi_iab
};

}  // namespace astdyn::math

#endif  // ASTDYN_MATH_TENSOR3_HPP
