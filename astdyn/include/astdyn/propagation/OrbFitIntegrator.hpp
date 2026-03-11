#ifndef ASTDYN_ORBFIT_INTEGRATOR_HPP
#define ASTDYN_ORBFIT_INTEGRATOR_HPP

#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/orbfit/integrator.h"
#include "astdyn/core/Types.hpp"
#include <Eigen/Dense>

namespace astdyn::propagation {

/**
 * @brief Adaptive Dormand-Prince RK45 integrator wrapped from OrbFit.
 * 
 * Uses orbfit's `propagate` internal logic wrapped into AstDyn's `Integrator` interface.
 */
class OrbFitDPIntegrator : public Integrator {
public:
    explicit OrbFitDPIntegrator(double initial_step = 0.1,
                                double tolerance = 1e-10,
                                double min_step = 1e-6,
                                double max_step = 5.0)
        : h_init_(initial_step), tol_(tolerance), h_min_(min_step), h_max_(max_step) {}

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

    std::vector<Eigen::VectorXd> integrate_at(const DerivativeFunction& f,
                                              const Eigen::VectorXd& y0,
                                              double t0,
                                              const std::vector<double>& t_targets) override;

private:
    double h_init_;
    double tol_;
    double h_min_;
    double h_max_;
};

/**
 * @brief Fixed-step RK4 integrator modeled after orbfit's `rk4_step`.
 */
class OrbFitRK4Integrator : public Integrator {
public:
    explicit OrbFitRK4Integrator(double step_size = 0.1) : h_(step_size) {}

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

    std::vector<Eigen::VectorXd> integrate_at(const DerivativeFunction& f,
                                              const Eigen::VectorXd& y0,
                                              double t0,
                                              const std::vector<double>& t_targets) override;

private:
    double h_;
};

} // namespace astdyn::propagation

#endif // ASTDYN_ORBFIT_INTEGRATOR_HPP
