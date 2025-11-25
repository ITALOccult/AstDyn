/**
 * @file LinearAlgebra.cpp
 * @brief Implementation of advanced linear algebra operations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 */

#include <astdyn/math/LinearAlgebra.hpp>
#include <astdyn/core/Constants.hpp>
#include <Eigen/SVD>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>

namespace astdyn {
namespace math {

// ============================================================================
// QR Decomposition
// ============================================================================

QRDecomposition qr_decompose(const MatrixXd& matrix) {
    Eigen::HouseholderQR<MatrixXd> qr(matrix);
    
    QRDecomposition result;
    result.Q = qr.householderQ();
    result.R = qr.matrixQR().triangularView<Eigen::Upper>();
    result.success = true;
    
    return result;
}

std::optional<VectorXd> qr_solve(const MatrixXd& matrix, const VectorXd& rhs) {
    if (matrix.rows() != rhs.size()) {
        return std::nullopt;
    }
    
    Eigen::ColPivHouseholderQR<MatrixXd> qr(matrix);
    if (!qr.isInvertible()) {
        return std::nullopt;
    }
    
    return qr.solve(rhs);
}

std::optional<MatrixXd> qr_invert(const MatrixXd& matrix) {
    if (matrix.rows() != matrix.cols()) {
        return std::nullopt;
    }
    
    Eigen::ColPivHouseholderQR<MatrixXd> qr(matrix);
    if (!qr.isInvertible()) {
        return std::nullopt;
    }
    
    MatrixXd identity = MatrixXd::Identity(matrix.rows(), matrix.cols());
    return qr.solve(identity);
}

// ============================================================================
// Singular Value Decomposition
// ============================================================================

SVDResult svd_decompose(const MatrixXd& matrix) {
    Eigen::JacobiSVD<MatrixXd> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    SVDResult result;
    result.U = svd.matrixU();
    result.singular_values = svd.singularValues();
    result.V = svd.matrixV();
    result.success = true;
    
    return result;
}

MatrixXd pseudoinverse(const MatrixXd& matrix, double tolerance) {
    Eigen::JacobiSVD<MatrixXd> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    VectorXd singular_values_inv = svd.singularValues();
    double threshold = tolerance * singular_values_inv(0);
    
    for (int i = 0; i < singular_values_inv.size(); ++i) {
        if (singular_values_inv(i) > threshold) {
            singular_values_inv(i) = 1.0 / singular_values_inv(i);
        } else {
            singular_values_inv(i) = 0.0;
        }
    }
    
    return svd.matrixV() * singular_values_inv.asDiagonal() * svd.matrixU().transpose();
}

// ============================================================================
// Eigenvalue Problems
// ============================================================================

EigenDecomposition eigen_symmetric(const MatrixXd& matrix) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(matrix);
    
    EigenDecomposition result;
    result.eigenvalues = solver.eigenvalues();
    result.eigenvectors = solver.eigenvectors();
    result.success = (solver.info() == Eigen::Success);
    
    return result;
}

EigenDecomposition eigen_general(const MatrixXd& matrix) {
    Eigen::EigenSolver<MatrixXd> solver(matrix);
    
    EigenDecomposition result;
    result.eigenvalues = solver.eigenvalues().real();
    result.eigenvectors = solver.eigenvectors().real();
    result.success = (solver.info() == Eigen::Success);
    
    return result;
}

std::optional<std::pair<double, VectorXd>> 
power_iteration(const MatrixXd& matrix, int max_iterations, double tolerance) {
    if (matrix.rows() != matrix.cols()) {
        return std::nullopt;
    }
    
    int n = matrix.rows();
    VectorXd v = VectorXd::Random(n);
    v.normalize();
    
    double eigenvalue = 0.0;
    double prev_eigenvalue;
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        VectorXd v_new = matrix * v;
        
        prev_eigenvalue = eigenvalue;
        eigenvalue = v.dot(v_new);
        
        v_new.normalize();
        v = v_new;
        
        if (std::abs(eigenvalue - prev_eigenvalue) < tolerance) {
            return std::make_pair(eigenvalue, v);
        }
    }
    
    return std::nullopt;  // Did not converge
}

// ============================================================================
// Least Squares
// ============================================================================

std::optional<VectorXd> least_squares(const MatrixXd& matrix, const VectorXd& rhs) {
    if (matrix.rows() != rhs.size()) {
        return std::nullopt;
    }
    
    Eigen::ColPivHouseholderQR<MatrixXd> qr(matrix);
    return qr.solve(rhs);
}

std::optional<VectorXd> weighted_least_squares(
    const MatrixXd& matrix, 
    const VectorXd& rhs,
    const MatrixXd& weight_matrix
) {
    // Transform to standard least squares: (W^{1/2} * A) * x = W^{1/2} * b
    Eigen::LLT<MatrixXd> llt(weight_matrix);
    if (llt.info() != Eigen::Success) {
        return std::nullopt;
    }
    
    MatrixXd sqrt_W = llt.matrixL();
    MatrixXd A_weighted = sqrt_W * matrix;
    VectorXd b_weighted = sqrt_W * rhs;
    
    return least_squares(A_weighted, b_weighted);
}

std::optional<std::pair<double, double>> 
linear_regression(const VectorXd& x, const VectorXd& y) {
    if (x.size() != y.size() || x.size() < 2) {
        return std::nullopt;
    }
    
    int n = x.size();
    MatrixXd A(n, 2);
    A.col(0) = VectorXd::Ones(n);
    A.col(1) = x;
    
    auto solution = least_squares(A, y);
    if (!solution) {
        return std::nullopt;
    }
    
    return std::make_pair((*solution)(1), (*solution)(0));  // (slope, intercept)
}

std::optional<std::tuple<double, double, double>>
quadratic_fit(const VectorXd& x, const VectorXd& y) {
    if (x.size() != y.size() || x.size() < 3) {
        return std::nullopt;
    }
    
    int n = x.size();
    MatrixXd A(n, 3);
    A.col(0) = VectorXd::Ones(n);
    A.col(1) = x;
    A.col(2) = x.array().square();
    
    auto solution = least_squares(A, y);
    if (!solution) {
        return std::nullopt;
    }
    
    return std::make_tuple((*solution)(0), (*solution)(1), (*solution)(2));
}

// ============================================================================
// LU Decomposition
// ============================================================================

LUDecomposition lu_decompose(const MatrixXd& matrix) {
    Eigen::FullPivLU<MatrixXd> lu(matrix);
    
    LUDecomposition result;
    result.L = lu.matrixLU().triangularView<Eigen::StrictlyLower>();
    result.L.diagonal() = VectorXd::Ones(matrix.rows());
    result.U = lu.matrixLU().triangularView<Eigen::Upper>();
    result.P = lu.permutationP();
    result.success = true;
    
    return result;
}

double determinant(const MatrixXd& matrix) {
    return matrix.determinant();
}

// ============================================================================
// Orthogonalization
// ============================================================================

MatrixXd gram_schmidt(const MatrixXd& vectors) {
    MatrixXd orthogonal = vectors;
    int n = vectors.cols();
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            double projection = orthogonal.col(i).dot(orthogonal.col(j));
            orthogonal.col(i) -= projection * orthogonal.col(j);
        }
        orthogonal.col(i).normalize();
    }
    
    return orthogonal;
}

MatrixXd modified_gram_schmidt(const MatrixXd& vectors) {
    MatrixXd orthogonal = vectors;
    int n = vectors.cols();
    
    for (int i = 0; i < n; ++i) {
        orthogonal.col(i).normalize();
        
        for (int j = i + 1; j < n; ++j) {
            double projection = orthogonal.col(j).dot(orthogonal.col(i));
            orthogonal.col(j) -= projection * orthogonal.col(i);
        }
    }
    
    return orthogonal;
}

// ============================================================================
// Matrix Exponential and Logarithm
// ============================================================================

MatrixXd matrix_exp(const MatrixXd& matrix) {
    return matrix.exp();
}

std::optional<MatrixXd> matrix_log(const MatrixXd& matrix) {
    // Check if matrix is positive definite
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(matrix);
    if (solver.info() != Eigen::Success) {
        return std::nullopt;
    }
    
    if ((solver.eigenvalues().array() <= 0).any()) {
        return std::nullopt;  // Not positive definite
    }
    
    return matrix.log();
}

// ============================================================================
// Special Matrix Operations
// ============================================================================

MatrixXd extract_block(const MatrixXd& matrix,
                       int start_row, int start_col,
                       int num_rows, int num_cols) {
    return matrix.block(start_row, start_col, num_rows, num_cols);
}

MatrixXd block_diagonal(const std::vector<MatrixXd>& blocks) {
    int total_rows = 0;
    int total_cols = 0;
    
    for (const auto& block : blocks) {
        total_rows += block.rows();
        total_cols += block.cols();
    }
    
    MatrixXd result = MatrixXd::Zero(total_rows, total_cols);
    
    int current_row = 0;
    int current_col = 0;
    
    for (const auto& block : blocks) {
        result.block(current_row, current_col, block.rows(), block.cols()) = block;
        current_row += block.rows();
        current_col += block.cols();
    }
    
    return result;
}

MatrixXd kronecker_product(const MatrixXd& A, const MatrixXd& B) {
    int m = A.rows();
    int n = A.cols();
    int p = B.rows();
    int q = B.cols();
    
    MatrixXd result(m * p, n * q);
    
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            result.block(i * p, j * q, p, q) = A(i, j) * B;
        }
    }
    
    return result;
}

} // namespace math
} // namespace astdyn
