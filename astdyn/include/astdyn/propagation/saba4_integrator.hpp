/**
 * @file saba4_integrator.hpp
 * @brief SABA4 (Symmetric Adaptive BAsed, order 4) symplectic integrator
 * @author Hybrid design: RKF78 adaptivity + Gauss symplectic structure
 */

#ifndef SABA4_INTEGRATOR_HPP
#define SABA4_INTEGRATOR_HPP

#include "astdyn/propagation/Integrator.hpp"
#include <Eigen/Dense>
#include <vector>

namespace astdyn::propagation {

class SABA4Integrator : public Integrator {
public:

    SABA4Integrator(double initial_step = 1.0,
                    double tolerance = 1e-6,
                    double min_step = 1e-6,
                    double max_step = 100.0);

    Eigen::VectorXd integrate(const DerivativeFunction& f,
                              const Eigen::VectorXd& y0,
                              double t0, double tf) override;

    void integrate_steps(const DerivativeFunction& f,
                         const Eigen::VectorXd& y0,
                         double t0, double tf,
                         std::vector<double>& t_out,
                         std::vector<Eigen::VectorXd>& y_out) override;

    // Use base class statistics() and reset_statistics()

private:
    static constexpr int num_stages_ = 7;
    static const double c_coeffs_[num_stages_];
    static const double d_coeffs_[num_stages_];

    Eigen::VectorXd saba4_step(const DerivativeFunction& f,
                               const Eigen::VectorXd& y,
                               double t, double h);

    Eigen::VectorXd saba2_step(const DerivativeFunction& f,
                               const Eigen::VectorXd& y,
                               double t, double h);

    void split_state(const Eigen::VectorXd& y,
                     Eigen::VectorXd& q, Eigen::VectorXd& p) const;
    Eigen::VectorXd join_state(const Eigen::VectorXd& q,
                               const Eigen::VectorXd& p) const;

    Eigen::VectorXd compute_force(const DerivativeFunction& f,
                                  double t,
                                  const Eigen::VectorXd& q,
                                  const Eigen::VectorXd& p) const;

    double adapt_step_size(double h_current,
                          double error_estimate,
                          bool last_accepted);

    double compute_energy(const Eigen::VectorXd& y) const;

    double h_initial_;
    double tolerance_;
    double h_min_;
    double h_max_;

    double initial_energy_ = 0.0;
};

} // namespace astdyn::propagation

#endif // SABA4_INTEGRATOR_HPP
