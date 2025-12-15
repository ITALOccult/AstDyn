/**
 * @file DifferentialCorrector.hpp
 * @brief Differential corrections for orbit determination
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 * 
 * Implements least squares differential corrections to refine orbital elements
 * from astrometric observations.
 * 
 * Algorithm:
 * 1. Compute residuals O-C for current orbit
 * 2. Compute design matrix A = ∂(O-C)/∂x₀ using STM
 * 3. Solve normal equations: (AᵀWA)Δx = AᵀWb
 * 4. Update orbit: x₀_new = x₀ + Δx
 * 5. Iterate until convergence
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP
#define ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include <memory>
#include <functional>

namespace astdyn::orbit_determination {

/**
 * @brief Settings for differential corrections
 */
struct DifferentialCorrectorSettings {
    int max_iterations = 20;             ///< Maximum iterations
    double convergence_tolerance = 1e-6; ///< Convergence threshold [AU]
    double outlier_sigma = 3.0;          ///< Target Sigma threshold (default)
    
    // Carpentry Settings (Iterative Rejection)
    double outlier_max_sigma = 10.0;     ///< Starting loose sigma
    double outlier_min_sigma = 3.0;      ///< Final tight sigma
    
    bool reject_outliers = true;         ///< Automatically reject outliers
    bool compute_covariance = true;      ///< Compute covariance matrix
    bool verbose = false;                ///< Print iteration details
};

/**
 * @brief Result of differential corrections
 */
struct DifferentialCorrectorResult {
    // Final orbit
    astdyn::propagation::CartesianElements final_state;
    bool converged;
    int iterations;
    
    // Quality metrics
    ResidualStatistics statistics;
    std::vector<ObservationResidual> residuals;
    
    // Uncertainty
    astdyn::Matrix6d covariance;                 ///< Covariance matrix [AU², AU²/day², ...]
    Eigen::VectorXd formal_uncertainties; ///< sqrt(diag(covariance))
    
    // Correlation matrix (normalized covariance)
    astdyn::Matrix6d correlation;
    
    // Normal matrix and its inverse
    astdyn::Matrix6d normal_matrix;              ///< AᵀWA
    astdyn::Matrix6d normal_matrix_inv;          ///< (AᵀWA)⁻¹
    
    // Convergence history
    std::vector<double> rms_history;     ///< RMS at each iteration
    std::vector<double> correction_norm; ///< ||Δx|| at each iteration
    
    /**
     * @brief Get 1-sigma uncertainty for parameter
     */
    double get_uncertainty(int param_index) const {
        return (param_index < 6) ? formal_uncertainties[param_index] : 0.0;
    }
    
    /**
     * @brief Print summary to stdout
     */
    void print_summary() const;
};

/**
 * @brief Differential corrector for orbit determination
 * 
 * Uses weighted least squares with state transition matrix to refine
 * an initial orbit estimate from astrometric observations.
 */
class DifferentialCorrector {
public:
    /**
     * @brief Constructor
     * 
     * @param residual_calc Residual calculator
     * @param stm_computer STM computer
     */
    DifferentialCorrector(
        std::shared_ptr<ResidualCalculator> residual_calc,
        std::shared_ptr<StateTransitionMatrix> stm_computer);
    
    /**
     * @brief Fit orbit to observations
     * 
     * @param observations Optical observations
     * @param initial_guess Initial orbital state
     * @param settings Iteration settings
     * @return Fit result with refined orbit and uncertainties
     */
    DifferentialCorrectorResult fit(
        const std::vector<astdyn::observations::OpticalObservation>& observations,
        const astdyn::propagation::CartesianElements& initial_guess,
        const DifferentialCorrectorSettings& settings = {});
    
    /**
     * @brief Perform single differential correction iteration
     * 
     * @param observations Observations
     * @param current_state Current orbital state
     * @param[out] correction State correction Δx
     * @param[out] residuals Computed residuals
     * @return True if iteration succeeded
     */
    bool iteration(
        const std::vector<astdyn::observations::OpticalObservation>& observations,
        const astdyn::propagation::CartesianElements& current_state,
        Eigen::VectorXd& correction,
        std::vector<ObservationResidual>& residuals);
    
    /**
     * @brief Compute covariance matrix from final fit
     * 
     * @param observations Observations used in fit
     * @param final_state Final orbital state
     * @param residuals Final residuals
     * @return Covariance matrix
     */
    astdyn::Matrix6d compute_covariance(
        const std::vector<astdyn::observations::OpticalObservation>& observations,
        const astdyn::propagation::CartesianElements& final_state,
        const std::vector<ObservationResidual>& residuals);
    
    /**
     * @brief Set callback for iteration progress
     * 
     * Callback signature: void(int iteration, const ResidualStatistics& stats)
     */
    using IterationCallback = std::function<void(int, const ResidualStatistics&)>;
    void set_iteration_callback(IterationCallback callback) {
        iteration_callback_ = callback;
    }

private:
    /**
     * @brief Build design matrix A (observation partials)
     * 
     * A is (2*N_obs x 6) matrix where each row pair represents:
     * ∂(RA,Dec)/∂x₀ = ∂(RA,Dec)/∂x * Φ(t,t₀)
     * 
     * @param observations Observations
     * @param state Current state
     * @return Design matrix and valid observation indices
     */
    struct DesignMatrixResult {
        Eigen::MatrixXd A;               ///< Design matrix (2N x 6)
        Eigen::VectorXd b;               ///< Residual vector (2N x 1)
        Eigen::VectorXd weights;         ///< Weight vector (2N x 1)
        std::vector<size_t> valid_indices; ///< Indices of valid observations
    };
    
    DesignMatrixResult build_design_matrix(
        const std::vector<astdyn::observations::OpticalObservation>& observations,
        const astdyn::propagation::CartesianElements& state,
        const std::vector<ObservationResidual>& residuals);
    
    /**
     * @brief Solve normal equations
     * 
     * (AᵀWA)Δx = AᵀWb
     * 
     * @param A Design matrix
     * @param b Residual vector
     * @param W Weight matrix (diagonal)
     * @param[out] normal_matrix Normal matrix AᵀWA
     * @param[out] normal_inv Inverse (AᵀWA)⁻¹
     * @return State correction Δx, or nullopt if singular
     */
    std::optional<Eigen::VectorXd> solve_normal_equations(
        const Eigen::MatrixXd& A,
        const Eigen::VectorXd& b,
        const Eigen::VectorXd& W,
        astdyn::Matrix6d& normal_matrix,
        astdyn::Matrix6d& normal_inv);
    
    /**
     * @brief Compute correlation matrix from covariance
     * 
     * corr(i,j) = cov(i,j) / sqrt(cov(i,i) * cov(j,j))
     */
    astdyn::Matrix6d compute_correlation(const astdyn::Matrix6d& covariance) const;
    
    /**
     * @brief Check convergence criteria
     */
    bool check_convergence(
        const Eigen::VectorXd& correction,
        double tolerance) const;

private:
    std::shared_ptr<ResidualCalculator> residual_calc_;
    std::shared_ptr<StateTransitionMatrix> stm_computer_;
    
    IterationCallback iteration_callback_;
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP
