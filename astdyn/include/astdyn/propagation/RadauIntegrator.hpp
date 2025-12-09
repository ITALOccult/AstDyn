/**
 * @file RadauIntegrator.hpp
 * @brief Radau IIA implicit integrator (15th order)
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * ⚠️  WARNING: NOT YET OPTIMIZED FOR PRODUCTION USE ⚠️
 * 
 * This integrator is currently SLOW and NOT RECOMMENDED for general use.
 * 
 * CURRENT STATUS:
 * - ❌ Too slow for non-stiff problems (100-1000× slower than RKF78)
 * - ❌ Newton solver not optimized
 * - ❌ No step size prediction
 * - ⚠️  Use ONLY for specific stiff problems
 * 
 * RECOMMENDED ALTERNATIVES:
 * - For general propagation: RKF78 (fast, accurate, adaptive)
 * - For orbit determination: RKF78 + STMPropagator
 * - For long-term: Gauss (symplectic, energy-conserving)
 * 
 * USE RADAU15 ONLY IF:
 * - Problem is genuinely stiff (rare in orbital mechanics)
 * - Maximum precision required (> 1e-13)
 * - You understand implicit methods and their cost
 * 
 * FUTURE WORK:
 * - Optimize Newton solver
 * - Add step size prediction
 * - Implement variable-order scheme
 * - Benchmark against DOPRI/RADAU from literature
 * 
 * References:
 * - Everhart, E. (1985) "An efficient integrator that uses Gauss-Radau spacings"
 * - Hairer & Wanner (1996) "Solving ODEs II: Stiff and DAE Problems"
 */

#ifndef ASTDYN_RADAU_INTEGRATOR_HPP
#define ASTDYN_RADAU_INTEGRATOR_HPP

#include "Integrator.hpp"
#include <Eigen/Dense>

namespace astdyn::propagation {

/**
 * @brief Radau IIA implicit integrator (15th order)
 * 
 * Radau IIA is a high-order implicit Runge-Kutta method with excellent
 * stability properties, particularly suited for stiff problems and
 * long-term orbit propagation.
 * 
 * Key features:
 * - Order 15 (very high accuracy)
 * - A-stable (excellent for stiff problems)
 * - Implicit (requires solving nonlinear system at each step)
 * - Adaptive step size control
 * 
 * The method uses 8 stages with Radau quadrature points.
 * At each step, it solves a nonlinear system using simplified Newton iteration.
 * 
 * References:
 * - Hairer & Wanner (1996) "Solving Ordinary Differential Equations II"
 * - Everhart (1985) "An efficient integrator that uses Gauss-Radau spacings"
 */
class RadauIntegrator : public Integrator {
public:
    /**
     * @brief Construct Radau15 integrator
     * 
     * @param initial_step Initial step size guess
     * @param tolerance Relative error tolerance (default 1e-13)
     * @param min_step Minimum allowed step size (default 1e-8)
     * @param max_step Maximum allowed step size (default 100.0)
     * @param max_newton_iter Maximum Newton iterations per step (default 4, optimized)
     */
    explicit RadauIntegrator(double initial_step,
                            double tolerance = 1e-13,
                            double min_step = 1e-8,
                            double max_step = 100.0,
                            int max_newton_iter = 7);
    
    Eigen::VectorXd integrate(const DerivativeFunction& f,
                             const Eigen::VectorXd& y0,
                             double t0,
                             double tf) override;
    
    void integrate_steps(const DerivativeFunction& f,
                        const Eigen::VectorXd& y0,
                        double t0,
                        double tf,
                        std::vector<double>& t_out,
                        std::vector<Eigen::VectorXd>& y_out) override;
    
    /**
     * @brief Adaptive Radau step with error control
     * 
     * @param f Derivative function
     * @param jac Jacobian function df/dy (optional, computed numerically if nullptr)
     * @param t Current time [in/out]
     * @param y Current state [in/out]
     * @param h Current step size [in/out]
     * @param t_target Target time
     * @return true if step accepted
     */
    bool adaptive_step(const DerivativeFunction& f,
                      std::function<Eigen::MatrixXd(double, const Eigen::VectorXd&)> jac,
                      double& t,
                      Eigen::VectorXd& y,
                      double& h,
                      double t_target);
    
    void set_tolerance(double tol) { tolerance_ = tol; }
    double tolerance() const { return tolerance_; }
    
    void set_max_newton_iter(int max_iter) { max_newton_iter_ = max_iter; }
    int max_newton_iter() const { return max_newton_iter_; }
    
private:
    double h_initial_;      ///< Initial step size
    double tolerance_;      ///< Relative error tolerance
    double h_min_;          ///< Minimum step size
    double h_max_;          ///< Maximum step size
    int max_newton_iter_;   ///< Max Newton iterations
    
    // Radau IIA coefficients (8 stages, order 15)
    static constexpr int num_stages_ = 8;
    static const double c_[num_stages_];      ///< Time nodes (Radau points)
    static const double a_[num_stages_][num_stages_];  ///< Runge-Kutta matrix
    static const double b_[num_stages_];      ///< Weights
    static const double b_hat_[num_stages_];  ///< Error estimator weights
    
    /**
     * @brief Compute numerical Jacobian
     */
    Eigen::MatrixXd numerical_jacobian(const DerivativeFunction& f,
                                       double t,
                                       const Eigen::VectorXd& y);
    
    /**
     * @brief Solve implicit Radau system using simplified Newton
     */
    bool solve_implicit_system(const DerivativeFunction& f,
                               const Eigen::MatrixXd& jacobian,
                               double t,
                               const Eigen::VectorXd& y,
                               double h,
                               std::vector<Eigen::VectorXd>& k);
};

} // namespace astdyn::propagation

#endif // ASTDYN_RADAU_INTEGRATOR_HPP
