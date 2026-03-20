/**
 * @file LeastSquaresFitter.hpp
 * @brief Least squares orbit fitter for orbit determination
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * Implements differential correction algorithm:
 * 1. Propagate orbit with STM
 * 2. Compute residuals O-C
 * 3. Build design matrix A = ∂ρ/∂x₀
 * 4. Solve normal equations: (AᵀWA)δx = AᵀWΔρ
 * 5. Update elements: x₀ ← x₀ + δx
 * 6. Iterate until convergence
 */

#ifndef ASTDYN_LEAST_SQUARES_FITTER_HPP
#define ASTDYN_LEAST_SQUARES_FITTER_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <functional>
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/time/epoch.hpp"

namespace astdyn::orbit_determination {

/**
 * @brief Least squares fit result
 */
struct FitResult {
    physics::CartesianStateTyped<core::GCRF> state;     ///< Fitted state [r, v]
    astdyn::Matrix6d covariance;                        ///< Covariance matrix (SI units: m^2, m^2/s, etc.)
    std::vector<ObservationResidual> residuals;
    
    int num_iterations;
    int num_observations;
    int num_rejected;
    
    astrometry::Angle rms_ra;
    astrometry::Angle rms_dec;
    astrometry::Angle rms_total;
    
    double chi_squared;
    bool converged;
};

/**
 * @brief Least squares orbit fitter
 */
class LeastSquaresFitter {
public:
    /**
     * @brief Residual function signature
     */
    using ResidualFunction = std::function<std::vector<ObservationResidual>(
        const physics::CartesianStateTyped<core::GCRF>&, time::EpochTDB)>;
    
    /**
     * @brief STM function signature
     */
    using STMFunction = std::function<std::pair<
        physics::CartesianStateTyped<core::GCRF>,
        astdyn::Matrix6d
    >(const physics::CartesianStateTyped<core::GCRF>&, time::EpochTDB, time::EpochTDB)>;
    
    /**
     * @brief Construct fitter
     */
    LeastSquaresFitter();
    
    /**
     * @brief Fit orbit to observations
     */
    FitResult fit(
        const physics::CartesianStateTyped<core::GCRF>& initial_state,
        time::EpochTDB epoch,
        ResidualFunction residual_func,
        STMFunction stm_func
    );
    
    /**
     * @brief Set convergence tolerance
     */
    void set_tolerance(double tol) { tolerance_ = tol; }
    
    /**
     * @brief Set maximum iterations
     */
    void set_max_iterations(int max_iter) { max_iterations_ = max_iter; }
    
    /**
     * @brief Set outlier rejection threshold (sigma)
     */
    void set_outlier_threshold(double sigma) { outlier_threshold_ = sigma; }
    
    /**
     * @brief Enable/disable outlier rejection
     */
    void set_outlier_rejection(bool enable) { outlier_rejection_ = enable; }
    
private:
    double tolerance_ = 1e-6;
    int max_iterations_ = 10;
    double outlier_threshold_ = 3.0;  // 3-sigma
    bool outlier_rejection_ = true;
    
    /**
     * @brief Build design matrix A = ∂ρ/∂x₀
     * 
     * For each observation:
     *   A_row = (∂ρ/∂x) × Φ(t, t₀)
     * 
     * where ∂ρ/∂x is computed numerically from residual function
     */
    Eigen::MatrixXd build_design_matrix(
        const std::vector<ObservationResidual>& residuals,
        const physics::CartesianStateTyped<core::GCRF>& state,
        time::EpochTDB epoch,
        STMFunction stm_func
    );
    
    Eigen::Vector<double, 6> solve_normal_equations(
        const Eigen::MatrixXd& A,
        const Eigen::VectorXd& residuals,
        const Eigen::VectorXd& weights,
        astdyn::Matrix6d& covariance
    );
    
    /**
     * @brief Reject outliers (3-sigma)
     */
    int reject_outliers(std::vector<ObservationResidual>& residuals);
    
    /**
     * @brief Compute RMS
     */
    void compute_statistics(
        const std::vector<ObservationResidual>& residuals,
        FitResult& result
    );
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_LEAST_SQUARES_FITTER_HPP
