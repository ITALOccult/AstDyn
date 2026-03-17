/**
 * @file AASIntegrator.hpp
 * @brief AAS (AstDyn-Adaptive Symplectic) integrator
 * 
 * Adaptive symplectic integrator of order 4 (Yoshida-based)
 * with a step size metric based on local potential curvature.
 */

#ifndef ASTDYN_AAS_INTEGRATOR_HPP
#define ASTDYN_AAS_INTEGRATOR_HPP

#include "astdyn/propagation/Integrator.hpp"
#include <Eigen/Dense>
#include <vector>

namespace astdyn::propagation {

/**
 * @brief AstDyn-Adaptive Symplectic (AAS) Integrator
 * 
 * An innovative symplectic integrator that adapts the step size based on
 * the local Hessian of the potential (spatial curvature) rather than
 * truncation error, preserving symplectic invariants while optimizing speed.
 */
class AASIntegrator : public Integrator {
public:
    /**
     * @brief Construct AAS Integrator
     * 
     * @param precision Target precision metric (influences ds)
     * @param mu Gravitational parameter of the central body [AU^3/day^2]
     * @param J2 Zonal harmonic coefficient (optional)
     * @param R_eq Equatorial radius for J2 calculation [AU]
     */
    AASIntegrator(double precision = 1e-4, 
                  double mu = 0.0002959122082855911,  // constants::GMS
                  double J2 = 0.0,
                  double R_eq = 0.004650467261);      // constants::R_SUN_AU

    Eigen::VectorXd integrate(const DerivativeFunction& f,
                              const Eigen::VectorXd& y0,
                              double t0, double tf) override;

    void integrate_steps(const DerivativeFunction& f,
                         const Eigen::VectorXd& y0,
                         double t0, double tf,
                         std::vector<double>& t_out,
                         std::vector<Eigen::VectorXd>& y_out) override;

    std::vector<Eigen::VectorXd> integrate_at(const DerivativeFunction& f,
                                             const Eigen::VectorXd& y0,
                                             double t0,
                                             const std::vector<double>& t_targets) override;

    void set_precision(double p) { precision_ = p; }
    void set_central_body(double mu, double J2 = 0.0, double R_eq = 6378137.0);

private:
    // Yoshida Order 4 Coefficients (Calculated in constructor)
    double w1, w0;
    double d1, d2;
    double c1, c2;

    double compute_modified_hamiltonian(const Eigen::VectorXd& q, const Eigen::VectorXd& p, double h) const;
    void symplectic_step(const DerivativeFunction& f, double t, Eigen::VectorXd& q, Eigen::VectorXd& p, Eigen::MatrixXd& phi, double ds);
    Eigen::Matrix3d compute_hessian(const Eigen::VectorXd& q) const;
    double compute_force_gradient(const Eigen::VectorXd& q) const;
    double compute_total_energy(const Eigen::VectorXd& q, const Eigen::VectorXd& p) const;
    double estimate_step_size(const Eigen::VectorXd& q, const Eigen::VectorXd& p, double target_dt) const;
    void update_phi_kick(const DerivativeFunction& f, double t, double step, const Eigen::VectorXd& q, const Eigen::VectorXd& p, Eigen::MatrixXd& phi);
    void split_state(const Eigen::VectorXd& y, Eigen::VectorXd& q, Eigen::VectorXd& p) const;
    Eigen::VectorXd join_state(const Eigen::VectorXd& q, const Eigen::VectorXd& p) const;
    Eigen::VectorXd finalize_state_phi(const Eigen::VectorXd& y0, const Eigen::VectorXd& q, const Eigen::VectorXd& p, const Eigen::MatrixXd& phi);

    double precision_;
    double mu_;
    double j2_;
    double r_eq_;
};

} // namespace astdyn::propagation

#endif // ASTDYN_AAS_INTEGRATOR_HPP
