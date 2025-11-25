/**
 * @file Integrator.hpp
 * @brief Numerical integrators for orbital propagation
 * 
 * This module provides numerical integration methods for solving
 * ordinary differential equations (ODEs) of the form:
 *   dy/dt = f(t, y)
 * 
 * Integrators implemented:
 * - RK4: Classic 4th-order Runge-Kutta (fixed step)
 * - RKF78: Runge-Kutta-Fehlberg 7(8) (adaptive step)
 */

#ifndef ASTDYN_INTEGRATOR_HPP
#define ASTDYN_INTEGRATOR_HPP

#include "astdyn/core/Types.hpp"
#include <Eigen/Dense>
#include <functional>
#include <vector>

namespace astdyn::propagation {

/**
 * @brief State derivative function signature
 * 
 * @param t Current time
 * @param y Current state vector
 * @return dy/dt State derivative
 */
using DerivativeFunction = std::function<Eigen::VectorXd(double t, const Eigen::VectorXd& y)>;

/**
 * @brief Integration statistics and diagnostics
 */
struct IntegrationStatistics {
    int num_steps = 0;           ///< Total steps taken
    int num_function_evals = 0;  ///< Total function evaluations
    int num_rejected_steps = 0;  ///< Steps rejected (adaptive only)
    double min_step_size = 0.0;  ///< Minimum step size used
    double max_step_size = 0.0;  ///< Maximum step size used
    double final_time = 0.0;     ///< Actual final time reached
    
    void reset() {
        num_steps = 0;
        num_function_evals = 0;
        num_rejected_steps = 0;
        min_step_size = 0.0;
        max_step_size = 0.0;
        final_time = 0.0;
    }
};

/**
 * @brief Base class for numerical integrators
 */
class Integrator {
public:
    virtual ~Integrator() = default;
    
    /**
     * @brief Integrate from t0 to tf
     * 
     * @param f Derivative function dy/dt = f(t, y)
     * @param y0 Initial state at t0
     * @param t0 Initial time
     * @param tf Final time
     * @return Final state at tf
     */
    virtual Eigen::VectorXd integrate(const DerivativeFunction& f,
                                      const Eigen::VectorXd& y0,
                                      double t0,
                                      double tf) = 0;
    
    /**
     * @brief Integrate and store intermediate steps
     * 
     * @param f Derivative function
     * @param y0 Initial state
     * @param t0 Initial time
     * @param tf Final time
     * @param t_out Output times
     * @param y_out Output states (resized to match t_out)
     */
    virtual void integrate_steps(const DerivativeFunction& f,
                                 const Eigen::VectorXd& y0,
                                 double t0,
                                 double tf,
                                 std::vector<double>& t_out,
                                 std::vector<Eigen::VectorXd>& y_out) = 0;
    
    /**
     * @brief Get integration statistics
     */
    const IntegrationStatistics& statistics() const { return stats_; }
    
    /**
     * @brief Reset statistics
     */
    void reset_statistics() { stats_.reset(); }
    
protected:
    IntegrationStatistics stats_;
};

/**
 * @brief Classic 4th-order Runge-Kutta integrator (fixed step)
 * 
 * RK4 is a robust fixed-step method with O(h^4) local truncation error.
 * Simple and stable for most orbital mechanics problems.
 * 
 * The method computes:
 *   k1 = f(t, y)
 *   k2 = f(t + h/2, y + h*k1/2)
 *   k3 = f(t + h/2, y + h*k2/2)
 *   k4 = f(t + h, y + h*k3)
 *   y(t+h) = y(t) + h*(k1 + 2*k2 + 2*k3 + k4)/6
 */
class RK4Integrator : public Integrator {
public:
    /**
     * @brief Construct RK4 integrator
     * 
     * @param step_size Fixed integration step size [same units as time]
     */
    explicit RK4Integrator(double step_size);
    
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
     * @brief Single RK4 step
     * 
     * @param f Derivative function
     * @param t Current time
     * @param y Current state
     * @param h Step size
     * @return New state y(t+h)
     */
    Eigen::VectorXd step(const DerivativeFunction& f,
                        double t,
                        const Eigen::VectorXd& y,
                        double h);
    
    void set_step_size(double h) { h_ = h; }
    double step_size() const { return h_; }
    
private:
    double h_; ///< Fixed step size
};

/**
 * @brief Runge-Kutta-Fehlberg 7(8) adaptive integrator
 * 
 * RKF78 uses a pair of 7th and 8th order formulas to estimate
 * local truncation error and adapt step size automatically.
 * 
 * This is a high-accuracy method suitable for long-term orbit
 * propagation where efficiency and accuracy are critical.
 * 
 * The method uses 13 function evaluations per step but achieves
 * very high accuracy with adaptive step control.
 * 
 * Reference: Fehlberg (1968) NASA TR R-287
 */
class RKF78Integrator : public Integrator {
public:
    /**
     * @brief Construct RKF78 integrator
    /**
     * @param initial_step Initial step size guess
     * @param tolerance Relative error tolerance (default 1e-12)
     * @param min_step Minimum allowed step size (default 1e-6 days, ~0.1 seconds)
     * @param max_step Maximum allowed step size (default 100.0 days)
     */
    explicit RKF78Integrator(double initial_step,
                            double tolerance = 1e-12,
                            double min_step = 1e-6,
                            double max_step = 100.0);
    
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
     * @brief Adaptive RKF78 step with error control
     * 
     * @param f Derivative function
     * @param t Current time [in/out]
     * @param y Current state [in/out]
     * @param h Current step size [in/out]
     * @param t_target Target time (don't overshoot)
     * @return true if step accepted, false if rejected
     */
    bool adaptive_step(const DerivativeFunction& f,
                      double& t,
                      Eigen::VectorXd& y,
                      double& h,
                      double t_target);
    
    void set_tolerance(double tol) { tolerance_ = tol; }
    double tolerance() const { return tolerance_; }
    
    void set_min_step(double h_min) { h_min_ = h_min; }
    void set_max_step(double h_max) { h_max_ = h_max; }
    
private:
    double h_initial_; ///< Initial step size
    double tolerance_; ///< Relative error tolerance
    double h_min_;     ///< Minimum step size
    double h_max_;     ///< Maximum step size
    
    // RKF78 Butcher tableau coefficients
    static const double a_[13][12];  // Matrix A
    static const double b7_[13];     // 7th order weights
    static const double b8_[13];     // 8th order weights
    static const double c_[13];      // Time nodes
};

} // namespace astdyn::propagation

#endif // ASTDYN_INTEGRATOR_HPP
