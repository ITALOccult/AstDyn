/**
 * @file MathUtils.hpp
 * @brief Mathematical utility functions for OrbFit
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Conversion from Fortran math_lib.f90 module.
 * Provides linear algebra operations, matrix decompositions,
 * and various mathematical utilities for orbit determination.
 */

#ifndef ASTDYN_MATH_MATHUTILS_HPP
#define ASTDYN_MATH_MATHUTILS_HPP

#include <astdyn/core/Types.hpp>
#include <astdyn/core/Constants.hpp>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <vector>
#include <optional>

namespace astdyn {
namespace math {

// ============================================================================
// Matrix Inversion and Decomposition
// ============================================================================

/**
 * @brief Cholesky decomposition of a symmetric positive-definite matrix
 * 
 * Computes the Cholesky factorization L*L^T = A where L is lower triangular.
 * This is the C++ equivalent of Fortran's tchol routine.
 * 
 * @param matrix Input symmetric positive-definite matrix
 * @return Lower triangular Cholesky factor L, or nullopt if not positive-definite
 */
std::optional<MatrixXd> cholesky_decompose(const MatrixXd& matrix);

/**
 * @brief Invert a symmetric positive-definite matrix via Cholesky
 * 
 * Computes the inverse using Cholesky decomposition.
 * More stable and efficient than direct inversion for SPD matrices.
 * Equivalent to Fortran's tchinv routine.
 * 
 * @param matrix Input symmetric positive-definite matrix
 * @return Inverse matrix, or nullopt if singular
 */
std::optional<MatrixXd> cholesky_invert(const MatrixXd& matrix);

/**
 * @brief Invert a general square matrix
 * 
 * Uses LU decomposition for general matrices.
 * 
 * @param matrix Input square matrix
 * @return Inverse matrix, or nullopt if singular
 */
std::optional<MatrixXd> matrix_invert(const MatrixXd& matrix);

/**
 * @brief Invert a 2x2 matrix (specialized, fast version)
 * 
 * Direct formula for 2x2 matrix inversion.
 * Equivalent to Fortran's inv22 routine.
 * 
 * @param matrix 2x2 matrix
 * @return Inverse matrix, or nullopt if determinant is zero
 */
std::optional<Eigen::Matrix2d> 
invert_2x2(const Eigen::Matrix2d& matrix);

// ============================================================================
// Norms and Distance Metrics
// ============================================================================

/**
 * @brief Compute weighted RMS norm: sqrt(x^T * W * x / n)
 * 
 * Computes the root-mean-square norm with weight matrix W.
 * Equivalent to Fortran's snorm function.
 * 
 * @param vector Input vector x
 * @param weight_matrix Weight matrix W (must be symmetric positive-definite)
 * @return Weighted RMS norm
 */
double weighted_rms_norm(const VectorXd& vector, const MatrixXd& weight_matrix);

/**
 * @brief Compute standard RMS norm: sqrt(x^T * x / n)
 * 
 * Simple root-mean-square norm without weights.
 * Equivalent to Fortran's snormd function.
 * 
 * @param vector Input vector
 * @return RMS norm
 */
double rms_norm(const VectorXd& vector);

/**
 * @brief Compute Euclidean (L2) norm
 * 
 * @param vector Input vector
 * @return L2 norm
 */
inline double euclidean_norm(const VectorXd& vector) {
    return vector.norm();
}

/**
 * @brief Compute maximum (L-infinity) norm
 * 
 * @param vector Input vector
 * @return Maximum absolute value
 */
inline double max_norm(const VectorXd& vector) {
    return vector.lpNorm<Eigen::Infinity>();
}

// ============================================================================
// Bilinear Forms and Products
// ============================================================================

/**
 * @brief Compute bilinear form: x^T * A * y
 * 
 * Equivalent to Fortran's bilin function.
 * 
 * @param x First vector
 * @param matrix Matrix A
 * @param y Second vector
 * @return Scalar result of bilinear form
 */
double bilinear_form(const VectorXd& x, const MatrixXd& matrix, const VectorXd& y);

/**
 * @brief Compute quadratic form: x^T * A * x
 * 
 * @param vector Vector x
 * @param matrix Symmetric matrix A
 * @return Scalar result of quadratic form
 */
inline double quadratic_form(const VectorXd& vector, const MatrixXd& matrix) {
    return bilinear_form(vector, matrix, vector);
}

// ============================================================================
// Matrix Conversions and Transformations
// ============================================================================

/**
 * @brief Convert covariance matrix between different coordinate systems
 * 
 * Transforms covariance: C_new = J * C_old * J^T
 * where J is the Jacobian of the transformation.
 * Equivalent to Fortran's convertcov routine.
 * 
 * @param covariance Original covariance matrix
 * @param jacobian Jacobian of coordinate transformation
 * @return Transformed covariance matrix
 */
MatrixXd convert_covariance(const MatrixXd& covariance, const MatrixXd& jacobian);

/**
 * @brief Convert normal matrix (inverse covariance) between coordinate systems
 * 
 * Transforms normal matrix: N_new = J^(-T) * N_old * J^(-1)
 * Equivalent to Fortran's convertnor routine.
 * 
 * @param normal_matrix Original normal matrix
 * @param jacobian Jacobian of coordinate transformation
 * @return Transformed normal matrix, or nullopt if Jacobian is singular
 */
std::optional<MatrixXd> convert_normal_matrix(
    const MatrixXd& normal_matrix, 
    const MatrixXd& jacobian
);

/**
 * @brief Propagate normal matrix using state transition matrix
 * 
 * For NxN normal matrix propagation.
 * Equivalent to Fortran's norprs routine.
 * 
 * @param normal_matrix Normal matrix at initial time
 * @param state_transition State transition matrix Phi
 * @return Propagated normal matrix
 */
MatrixXd propagate_normal_matrix(
    const MatrixXd& normal_matrix,
    const MatrixXd& state_transition
);

// ============================================================================
// Numerical Utilities
// ============================================================================

/**
 * @brief Safe division with underflow/overflow protection
 * 
 * @param numerator Numerator
 * @param denominator Denominator
 * @param default_value Value to return if division is unsafe
 * @return Result of division, or default_value
 */
double safe_divide(double numerator, double denominator, 
                   double default_value = 0.0);

/**
 * @brief Check if a matrix is symmetric
 * 
 * @param matrix Input matrix
 * @param tolerance Tolerance for symmetry check
 * @return true if symmetric within tolerance
 */
bool is_symmetric(const MatrixXd& matrix, 
                  double tolerance = constants::DEFAULT_TOLERANCE);

/**
 * @brief Check if a matrix is positive definite
 * 
 * Uses eigenvalue decomposition.
 * 
 * @param matrix Input symmetric matrix
 * @return true if all eigenvalues are positive
 */
bool is_positive_definite(const MatrixXd& matrix);

/**
 * @brief Compute matrix rank
 * 
 * Uses singular value decomposition.
 * 
 * @param matrix Input matrix
 * @param tolerance Threshold for considering singular values as zero
 * @return Rank of matrix
 */
int matrix_rank(const MatrixXd& matrix, 
                double tolerance = constants::SMALL);

/**
 * @brief Compute condition number of a matrix
 * 
 * Ratio of largest to smallest singular value.
 * 
 * @param matrix Input matrix
 * @return Condition number (inf if singular)
 */
double condition_number(const MatrixXd& matrix);

// ============================================================================
// Statistical Functions
// ============================================================================

/**
 * @brief Compute mean of vector elements
 * 
 * @param vector Input vector
 * @return Mean value
 */
inline double mean(const VectorXd& vector) {
    return vector.mean();
}

/**
 * @brief Compute standard deviation of vector elements
 * 
 * @param vector Input vector
 * @param sample Use sample standard deviation (n-1) if true
 * @return Standard deviation
 */
double standard_deviation(const VectorXd& vector, bool sample = true);

/**
 * @brief Compute median of vector elements
 * 
 * @param vector Input vector
 * @return Median value
 */
double median(VectorXd vector);  // Note: takes copy for sorting

/**
 * @brief Compute covariance between two vectors
 * 
 * @param x First vector
 * @param y Second vector
 * @return Covariance value
 */
double covariance(const VectorXd& x, const VectorXd& y);

/**
 * @brief Compute correlation coefficient between two vectors
 * 
 * @param x First vector
 * @param y Second vector
 * @return Correlation coefficient [-1, 1]
 */
double correlation(const VectorXd& x, const VectorXd& y);

// ============================================================================
// Angle Utilities
// ============================================================================

/**
 * @brief Normalize angle to [0, 2π) range
 * 
 * @param angle Angle in radians
 * @return Normalized angle in [0, 2π)
 */
double normalize_angle_positive(double angle);

/**
 * @brief Normalize angle to [-π, π) range
 * 
 * @param angle Angle in radians
 * @return Normalized angle in [-π, π)
 */
double normalize_angle_symmetric(double angle);

/**
 * @brief Compute angular difference (shortest path)
 * 
 * Returns angle in [-π, π] representing shortest rotation from angle1 to angle2.
 * 
 * @param angle1 First angle in radians
 * @param angle2 Second angle in radians
 * @return Angular difference in [-π, π]
 */
double angular_difference(double angle1, double angle2);

} // namespace math
} // namespace astdyn

#endif // ASTDYN_MATH_MATHUTILS_HPP
