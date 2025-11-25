/**
 * @file LinearAlgebra.hpp
 * @brief Advanced linear algebra operations for OrbFit
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Additional linear algebra routines beyond basic MathUtils.
 * Includes QR decomposition, eigenvalue problems, least squares, etc.
 */

#ifndef ASTDYN_MATH_LINEARALGEBRA_HPP
#define ASTDYN_MATH_LINEARALGEBRA_HPP

#include <astdyn/core/Types.hpp>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <optional>

namespace astdyn {
namespace math {

// ============================================================================
// QR Decomposition and Related Operations
// ============================================================================

/**
 * @brief QR decomposition result
 */
struct QRDecomposition {
    MatrixXd Q;  ///< Orthogonal matrix
    MatrixXd R;  ///< Upper triangular matrix
    bool success; ///< Decomposition successful
};

/**
 * @brief Perform QR decomposition: A = Q*R
 * 
 * @param matrix Input matrix
 * @return QR decomposition result
 */
QRDecomposition qr_decompose(const MatrixXd& matrix);

/**
 * @brief Solve linear system using QR decomposition
 * 
 * Solves A*x = b using QR factorization.
 * 
 * @param matrix Coefficient matrix A
 * @param rhs Right-hand side vector b
 * @return Solution vector x, or nullopt if singular
 */
std::optional<VectorXd> qr_solve(const MatrixXd& matrix, const VectorXd& rhs);

/**
 * @brief Invert matrix using QR decomposition
 * 
 * @param matrix Input matrix
 * @return Inverse matrix, or nullopt if singular
 */
std::optional<MatrixXd> qr_invert(const MatrixXd& matrix);

// ============================================================================
// Singular Value Decomposition (SVD)
// ============================================================================

/**
 * @brief SVD decomposition result
 */
struct SVDResult {
    MatrixXd U;           ///< Left singular vectors
    VectorXd singular_values; ///< Singular values (sorted descending)
    MatrixXd V;           ///< Right singular vectors
    bool success;         ///< Decomposition successful
};

/**
 * @brief Perform Singular Value Decomposition: A = U*S*V^T
 * 
 * @param matrix Input matrix
 * @return SVD result
 */
SVDResult svd_decompose(const MatrixXd& matrix);

/**
 * @brief Compute Moore-Penrose pseudoinverse using SVD
 * 
 * @param matrix Input matrix
 * @param tolerance Threshold for considering singular values as zero
 * @return Pseudoinverse matrix
 */
MatrixXd pseudoinverse(const MatrixXd& matrix, 
                       double tolerance = 1e-10);

// ============================================================================
// Eigenvalue and Eigenvector Problems
// ============================================================================

/**
 * @brief Eigenvalue decomposition result
 */
struct EigenDecomposition {
    VectorXd eigenvalues;      ///< Eigenvalues
    MatrixXd eigenvectors;     ///< Eigenvectors (columns)
    bool success;              ///< Decomposition successful
};

/**
 * @brief Compute eigenvalues and eigenvectors of a symmetric matrix
 * 
 * More efficient than general eigenvalue solver for symmetric matrices.
 * 
 * @param matrix Input symmetric matrix
 * @return Eigenvalue decomposition
 */
EigenDecomposition eigen_symmetric(const MatrixXd& matrix);

/**
 * @brief Compute eigenvalues and eigenvectors of a general matrix
 * 
 * @param matrix Input matrix
 * @return Eigenvalue decomposition (may have complex eigenvalues)
 */
EigenDecomposition eigen_general(const MatrixXd& matrix);

/**
 * @brief Find largest eigenvalue and corresponding eigenvector (power method)
 * 
 * Useful for dominant eigenvalue problems.
 * 
 * @param matrix Input matrix
 * @param max_iterations Maximum iterations
 * @param tolerance Convergence tolerance
 * @return Pair of (eigenvalue, eigenvector), or nullopt if failed
 */
std::optional<std::pair<double, VectorXd>> 
power_iteration(const MatrixXd& matrix, 
                int max_iterations = 1000,
                double tolerance = 1.0e-10);

// ============================================================================
// Least Squares Problems
// ============================================================================

/**
 * @brief Solve linear least squares problem: min ||A*x - b||^2
 * 
 * Uses QR decomposition for stable solution.
 * 
 * @param matrix Coefficient matrix A (m x n, m >= n)
 * @param rhs Right-hand side b (m-vector)
 * @return Least squares solution x (n-vector), or nullopt if failed
 */
std::optional<VectorXd> least_squares(const MatrixXd& matrix, const VectorXd& rhs);

/**
 * @brief Weighted least squares: min ||(A*x - b)^T * W * (A*x - b)||
 * 
 * @param matrix Coefficient matrix A
 * @param rhs Right-hand side b
 * @param weight_matrix Weight matrix W (symmetric positive-definite)
 * @return Weighted least squares solution, or nullopt if failed
 */
std::optional<VectorXd> weighted_least_squares(
    const MatrixXd& matrix, 
    const VectorXd& rhs,
    const MatrixXd& weight_matrix
);

/**
 * @brief Linear regression: fit y = a + b*x
 * 
 * Equivalent to Fortran's linfi3 routine (for 1D case).
 * 
 * @param x Independent variable
 * @param y Dependent variable
 * @return Pair of (slope, intercept), or nullopt if failed
 */
std::optional<std::pair<double, double>> 
linear_regression(const VectorXd& x, const VectorXd& y);

/**
 * @brief Quadratic fit: y = a + b*x + c*x^2
 * 
 * Equivalent to Fortran's quadratic_fit routine.
 * 
 * @param x Independent variable
 * @param y Dependent variable
 * @return Triplet (a, b, c), or nullopt if failed
 */
std::optional<std::tuple<double, double, double>>
quadratic_fit(const VectorXd& x, const VectorXd& y);

// ============================================================================
// Matrix Factorizations
// ============================================================================

/**
 * @brief LU decomposition result
 */
struct LUDecomposition {
    MatrixXd L;           ///< Lower triangular
    MatrixXd U;           ///< Upper triangular
    Eigen::PermutationMatrix<Eigen::Dynamic> P;  ///< Permutation matrix
    bool success;
};

/**
 * @brief Perform LU decomposition with partial pivoting: P*A = L*U
 * 
 * @param matrix Input matrix
 * @return LU decomposition result
 */
LUDecomposition lu_decompose(const MatrixXd& matrix);

/**
 * @brief Compute matrix determinant via LU decomposition
 * 
 * @param matrix Input square matrix
 * @return Determinant value
 */
double determinant(const MatrixXd& matrix);

// ============================================================================
// Orthogonalization
// ============================================================================

/**
 * @brief Gram-Schmidt orthogonalization
 * 
 * Orthogonalizes a set of vectors (matrix columns).
 * 
 * @param vectors Input matrix (vectors as columns)
 * @return Orthogonalized vectors
 */
MatrixXd gram_schmidt(const MatrixXd& vectors);

/**
 * @brief Modified Gram-Schmidt (more numerically stable)
 * 
 * @param vectors Input matrix (vectors as columns)
 * @return Orthogonalized vectors
 */
MatrixXd modified_gram_schmidt(const MatrixXd& vectors);

// ============================================================================
// Matrix Exponential and Logarithm
// ============================================================================

/**
 * @brief Compute matrix exponential: exp(A)
 * 
 * Uses Padé approximation.
 * 
 * @param matrix Input matrix
 * @return Matrix exponential
 */
MatrixXd matrix_exp(const MatrixXd& matrix);

/**
 * @brief Compute matrix logarithm: log(A)
 * 
 * Requires A to be positive definite.
 * 
 * @param matrix Input matrix (must be positive definite)
 * @return Matrix logarithm, or nullopt if not possible
 */
std::optional<MatrixXd> matrix_log(const MatrixXd& matrix);

// ============================================================================
// Special Matrix Operations
// ============================================================================

/**
 * @brief Extract block from matrix
 * 
 * @param matrix Input matrix
 * @param start_row Starting row index
 * @param start_col Starting column index
 * @param num_rows Number of rows in block
 * @param num_cols Number of columns in block
 * @return Block matrix
 */
MatrixXd extract_block(const MatrixXd& matrix,
                       int start_row, int start_col,
                       int num_rows, int num_cols);

/**
 * @brief Build block diagonal matrix from list of matrices
 * 
 * @param blocks Vector of matrices to place on diagonal
 * @return Block diagonal matrix
 */
MatrixXd block_diagonal(const std::vector<MatrixXd>& blocks);

/**
 * @brief Kronecker product of two matrices
 * 
 * @param A First matrix
 * @param B Second matrix
 * @return Kronecker product A ⊗ B
 */
MatrixXd kronecker_product(const MatrixXd& A, const MatrixXd& B);

} // namespace math
} // namespace astdyn

#endif // ASTDYN_MATH_LINEARALGEBRA_HPP
