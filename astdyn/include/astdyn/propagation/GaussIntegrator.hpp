/**
 * @file GaussIntegrator.hpp
 * @brief Gauss-Legendre implicit symplectic integrator
 * @author AstDyn Team
 * @date 2025-12-09
 */

#ifndef ASTDYN_GAUSS_INTEGRATOR_HPP
#define ASTDYN_GAUSS_INTEGRATOR_HPP

#include "Integrator.hpp"
#include <Eigen/Dense>

namespace astdyn::propagation {

/**
 * @brief Gauss-Legendre implicit symplectic integrator
 * 
 * Gauss-Legendre methods are implicit Runge-Kutta methods with
 * Gauss quadrature points. They are symplectic, meaning they
 * preserve the Hamiltonian structure of the system.
 * 
 * Key features:
 * - Symplectic (preserves energy in Hamiltonian systems)
 * - Order 2s for s stages
 * - Excellent for long-term integrations (no secular drift)
 * - Requires solving implicit system at each step
 * 
 * This implementation uses 4 stages (order 8).
 * 
 * References:
 * - Hairer, Lubich, Wanner (2006) "Geometric Numerical Integration"
 * - Sanz-Serna & Calvo (1994) "Numerical Hamiltonian Problems"
 */
class GaussIntegrator : public Integrator {
public:
    /**
     * @brief Construct Gauss integrator
     * 
     * @param initial_step Initial step size
     * @param tolerance Error tolerance (default 1e-12)
     * @param min_step Minimum step size (default 1e-8)
     * @param max_step Maximum step size (default 100.0)
     * @param max_newton_iter Maximum Newton iterations (default 10)
     */
    explicit GaussIntegrator(double initial_step,
                            double tolerance = 1e-12,
                            double min_step = 1e-8,
                            double max_step = 100.0,
                            int max_newton_iter = 10);
    
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
    
    void set_tolerance(double tol) { tolerance_ = tol; }
    double tolerance() const { return tolerance_; }
    
private:
    double h_initial_;
    double tolerance_;
    double h_min_;
    double h_max_;
    int max_newton_iter_;
    
    // Gauss-Legendre coefficients (4 stages, order 8)
    static constexpr int num_stages_ = 4;
    static const double c_[num_stages_];      ///< Gauss points
    static const double a_[num_stages_][num_stages_];  ///< RK matrix
    static const double b_[num_stages_];      ///< Weights
    
    /**
     * @brief Solve implicit Gauss system (simplified Newton)
     */
    bool solve_implicit_system(const DerivativeFunction& f,
                               double t,
                               const Eigen::VectorXd& y,
                               double h,
                               std::vector<Eigen::VectorXd>& k);
    
    /**
     * @brief Compute energy for Hamiltonian systems (monitoring)
     */
    double compute_energy(const Eigen::VectorXd& y) const;
};

} // namespace astdyn::propagation

#endif // ASTDYN_GAUSS_INTEGRATOR_HPP
