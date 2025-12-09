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

#include <Eigen/Dense>
#include <vector>
#include <functional>

namespace astdyn::orbit_determination {

/**
 * @brief Observation with residual
 */
struct ObservationResidual {
    double epoch_mjd;
    double ra_obs_deg;          ///< Observed RA [deg]
    double dec_obs_deg;         ///< Observed Dec [deg]
    double ra_comp_deg;         ///< Computed RA [deg]
    double dec_comp_deg;        ///< Computed Dec [deg]
    double ra_residual_arcsec;  ///< RA*(O-C) [arcsec]
    double dec_residual_arcsec; ///< Dec(O-C) [arcsec]
    double weight_ra;           ///< Weight (1/σ²)
    double weight_dec;
    bool rejected;              ///< Outlier flag
};

/**
 * @brief Least squares fit result
 */
struct FitResult {
    Eigen::Vector<double, 6> state;     ///< Fitted state [r, v]
    Eigen::Matrix<double, 6, 6> covariance; ///< Covariance matrix
    std::vector<ObservationResidual> residuals;
    
    int num_iterations;
    int num_observations;
    int num_rejected;
    
    double rms_ra_arcsec;
    double rms_dec_arcsec;
    double rms_total_arcsec;
    
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
     * 
     * Given state at epoch, returns vector of residuals
     * Format: [ra1, dec1, ra2, dec2, ...] in arcsec
     */
    using ResidualFunction = std::function<std::vector<ObservationResidual>(
        const Eigen::Vector<double, 6>&, double)>;
    
    /**
     * @brief STM function signature
     * 
     * Returns state and STM at target epoch
     */
    using STMFunction = std::function<std::pair<
        Eigen::Vector<double, 6>,
        Eigen::Matrix<double, 6, 6>
    >(const Eigen::Vector<double, 6>&, double, double)>;
    
    /**
     * @brief Construct fitter
     */
    LeastSquaresFitter();
    
    /**
     * @brief Fit orbit to observations
     * 
     * @param initial_state Initial guess [r, v]
     * @param epoch_mjd Epoch of initial state
     * @param residual_func Function to compute residuals
     * @param stm_func Function to propagate with STM
     * @return Fit result
     */
    FitResult fit(
        const Eigen::Vector<double, 6>& initial_state,
        double epoch_mjd,
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
        const Eigen::Vector<double, 6>& state,
        double epoch_mjd,
        STMFunction stm_func
    );
    
    /**
     * @brief Solve normal equations
     * 
     * (AᵀWA)δx = AᵀWΔρ
     */
    Eigen::Vector<double, 6> solve_normal_equations(
        const Eigen::MatrixXd& A,
        const Eigen::VectorXd& residuals,
        const Eigen::VectorXd& weights,
        Eigen::Matrix<double, 6, 6>& covariance
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
