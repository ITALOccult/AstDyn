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
    struct Statistics {
        int num_steps = 0;
        int num_rejected = 0;
        int num_accepted = 0;
        int num_function_evals = 0;
        double energy_drift = 0.0;
        double final_time = 0.0;
    };

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

    Statistics stats() const { return stats_; }
    void reset_stats() { stats_ = Statistics(); }

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
                     Eigen::Vector3d& q, Eigen::Vector3d& p) const;
    Eigen::VectorXd join_state(const Eigen::Vector3d& q,
                               const Eigen::Vector3d& p) const;

    Eigen::Vector3d compute_force(const DerivativeFunction& f,
                                  double t,
                                  const Eigen::Vector3d& q,
                                  const Eigen::Vector3d& p) const;

    double adapt_step_size(double h_current,
                          double error_estimate,
                          bool last_accepted);

    double compute_energy(const Eigen::VectorXd& y) const;

    double h_initial_;
    double tolerance_;
    double h_min_;
    double h_max_;

    Statistics stats_;
    double initial_energy_ = 0.0;
};

} // namespace astdyn::propagation

#endif // SABA4_INTEGRATOR_HPP
