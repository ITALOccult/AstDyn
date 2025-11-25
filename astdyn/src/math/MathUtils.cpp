/**
 * @file MathUtils.cpp
 * @brief Implementation of mathematical utility functions
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 */

#include <astdyn/math/MathUtils.hpp>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace astdyn {
namespace math {

// ============================================================================
// Matrix Inversion and Decomposition
// ============================================================================

std::optional<MatrixXd> cholesky_decompose(const MatrixXd& matrix) {
    if (matrix.rows() != matrix.cols()) {
        return std::nullopt;
    }
    
    Eigen::LLT<MatrixXd> llt(matrix);
    if (llt.info() != Eigen::Success) {
        return std::nullopt;
    }
    
    return llt.matrixL();
}

std::optional<MatrixXd> cholesky_invert(const MatrixXd& matrix) {
    if (matrix.rows() != matrix.cols()) {
        return std::nullopt;
    }
    
    Eigen::LLT<MatrixXd> llt(matrix);
    if (llt.info() != Eigen::Success) {
        return std::nullopt;
    }
    
    MatrixXd identity = MatrixXd::Identity(matrix.rows(), matrix.cols());
    return llt.solve(identity);
}

std::optional<MatrixXd> matrix_invert(const MatrixXd& matrix) {
    if (matrix.rows() != matrix.cols()) {
        return std::nullopt;
    }
    
    Eigen::FullPivLU<MatrixXd> lu(matrix);
    if (!lu.isInvertible()) {
        return std::nullopt;
    }
    
    return lu.inverse();
}

std::optional<Eigen::Matrix2d> 
invert_2x2(const Eigen::Matrix2d& matrix) {
    double det = matrix.determinant();
    
    if (std::abs(det) < 1e-15) {
        return std::nullopt;
    }
    
    Eigen::Matrix2d inverse;
    inverse(0, 0) =  matrix(1, 1) / det;
    inverse(0, 1) = -matrix(0, 1) / det;
    inverse(1, 0) = -matrix(1, 0) / det;
    inverse(1, 1) =  matrix(0, 0) / det;
    
    return inverse;
}

// ============================================================================
// Norms and Distance Metrics
// ============================================================================

double weighted_rms_norm(const VectorXd& vector, const MatrixXd& weight_matrix) {
    if (vector.size() == 0) {
        return 0.0;
    }
    
    if (weight_matrix.rows() != vector.size() || 
        weight_matrix.cols() != vector.size()) {
        throw std::invalid_argument("Weight matrix dimensions must match vector size");
    }
    
    double quadratic = vector.transpose() * weight_matrix * vector;
    return std::sqrt(quadratic / vector.size());
}

double rms_norm(const VectorXd& vector) {
    if (vector.size() == 0) {
        return 0.0;
    }
    
    return std::sqrt(vector.squaredNorm() / vector.size());
}

// ============================================================================
// Bilinear Forms and Products
// ============================================================================

double bilinear_form(const VectorXd& x, const MatrixXd& matrix, const VectorXd& y) {
    if (x.size() != matrix.rows() || y.size() != matrix.cols()) {
        throw std::invalid_argument("Incompatible dimensions for bilinear form");
    }
    
    return x.transpose() * matrix * y;
}

// ============================================================================
// Matrix Conversions and Transformations
// ============================================================================

MatrixXd convert_covariance(const MatrixXd& covariance, const MatrixXd& jacobian) {
    // C_new = J * C_old * J^T
    return jacobian * covariance * jacobian.transpose();
}

std::optional<MatrixXd> convert_normal_matrix(
    const MatrixXd& normal_matrix, 
    const MatrixXd& jacobian
) {
    // N_new = J^(-T) * N_old * J^(-1)
    // = (J^(-1))^T * N_old * J^(-1)
    
    auto jacobian_inv = matrix_invert(jacobian);
    if (!jacobian_inv) {
        return std::nullopt;
    }
    
    return jacobian_inv->transpose() * normal_matrix * (*jacobian_inv);
}

MatrixXd propagate_normal_matrix(
    const MatrixXd& normal_matrix,
    const MatrixXd& state_transition
) {
    // N(t) = Phi^(-T) * N(t0) * Phi^(-1)
    // where Phi is the state transition matrix
    
    auto result = convert_normal_matrix(normal_matrix, state_transition);
    if (!result) {
        throw std::runtime_error("State transition matrix is singular");
    }
    
    return *result;
}

// ============================================================================
// Numerical Utilities
// ============================================================================

double safe_divide(double numerator, double denominator, double default_value) {
    if (std::abs(denominator) < constants::EPSILON) {
        return default_value;
    }
    
    double result = numerator / denominator;
    
    // Check for overflow/underflow
    if (!std::isfinite(result)) {
        return default_value;
    }
    
    return result;
}

bool is_symmetric(const MatrixXd& matrix, double tolerance) {
    if (matrix.rows() != matrix.cols()) {
        return false;
    }
    
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = i + 1; j < matrix.cols(); ++j) {
            if (std::abs(matrix(i, j) - matrix(j, i)) > tolerance) {
                return false;
            }
        }
    }
    
    return true;
}

bool is_positive_definite(const MatrixXd& matrix) {
    if (!is_symmetric(matrix)) {
        return false;
    }
    
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(matrix);
    if (es.info() != Eigen::Success) {
        return false;
    }
    
    return (es.eigenvalues().array() > 0).all();
}

int matrix_rank(const MatrixXd& matrix, double tolerance) {
    Eigen::JacobiSVD<MatrixXd> svd(matrix);
    
    int rank = 0;
    double threshold = tolerance * svd.singularValues()(0);
    
    for (int i = 0; i < svd.singularValues().size(); ++i) {
        if (svd.singularValues()(i) > threshold) {
            rank++;
        }
    }
    
    return rank;
}

double condition_number(const MatrixXd& matrix) {
    Eigen::JacobiSVD<MatrixXd> svd(matrix);
    
    auto sv = svd.singularValues();
    if (sv.size() == 0 || sv(sv.size() - 1) < 1e-15) {
        return std::numeric_limits<double>::infinity();
    }
    
    return sv(0) / sv(sv.size() - 1);
}

// ============================================================================
// Statistical Functions
// ============================================================================

double standard_deviation(const VectorXd& vector, bool sample) {
    if (vector.size() <= 1) {
        return 0.0;
    }
    
    double mean_val = vector.mean();
    double variance = (vector.array() - mean_val).square().sum();
    
    int divisor = sample ? (vector.size() - 1) : vector.size();
    
    return std::sqrt(variance / divisor);
}

double median(VectorXd vector) {
    if (vector.size() == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    std::sort(vector.data(), vector.data() + vector.size());
    
    int mid = vector.size() / 2;
    if (vector.size() % 2 == 0) {
        return (vector(mid - 1) + vector(mid)) / 2.0;
    } else {
        return vector(mid);
    }
}

double covariance(const VectorXd& x, const VectorXd& y) {
    if (x.size() != y.size() || x.size() == 0) {
        throw std::invalid_argument("Vectors must have same non-zero size");
    }
    
    double mean_x = x.mean();
    double mean_y = y.mean();
    
    double cov = 0.0;
    for (int i = 0; i < x.size(); ++i) {
        cov += (x(i) - mean_x) * (y(i) - mean_y);
    }
    
    return cov / x.size();
}

double correlation(const VectorXd& x, const VectorXd& y) {
    double cov = covariance(x, y);
    double std_x = standard_deviation(x, false);
    double std_y = standard_deviation(y, false);
    
    if (std_x < 1e-15 || std_y < 1e-15) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    return cov / (std_x * std_y);
}

// ============================================================================
// Angle Utilities
// ============================================================================

double normalize_angle_positive(double angle) {
    double normalized = std::fmod(angle, constants::TWO_PI);
    if (normalized < 0.0) {
        normalized += constants::TWO_PI;
    }
    return normalized;
}

double normalize_angle_symmetric(double angle) {
    double normalized = std::fmod(angle + constants::PI, constants::TWO_PI);
    if (normalized < 0.0) {
        normalized += constants::TWO_PI;
    }
    return normalized - constants::PI;
}

double angular_difference(double angle1, double angle2) {
    double diff = angle2 - angle1;
    return normalize_angle_symmetric(diff);
}

} // namespace math
} // namespace astdyn
