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
     * @param mu Gravitational parameter of the central body [m^3/s^2]
     * @param J2 Zonal harmonic coefficient (optional, for Earth/Planets)
     * @param R_eq Equatorial radius for J2 calculation [m]
     */
    AASIntegrator(double precision = 1e-4, 
                  double mu = 1.32712440018e20,  // default: Sun
                  double J2 = 0.0,
                  double R_eq = 6378137.0);

    Eigen::VectorXd integrate(const DerivativeFunction& f,
                              const Eigen::VectorXd& y0,
                              double t0, double tf) override;

    void integrate_steps(const DerivativeFunction& f,
                         const Eigen::VectorXd& y0,
                         double t0, double tf,
                         std::vector<double>& t_out,
                         std::vector<Eigen::VectorXd>& y_out) override;

    void set_precision(double p) { precision_ = p; }
    void set_central_body(double mu, double J2 = 0.0, double R_eq = 6378137.0);

private:
    // Yoshida Order 4 Coefficients (Calculated in constructor)
    double w1, w0;
    double d1, d2;
    double c1, c2;

    double compute_modified_hamiltonian(const Eigen::VectorXd& q,
                                        const Eigen::VectorXd& p,
                                        double h) const;

    void symplectic_step(const DerivativeFunction& f,
                         double t,
                         Eigen::VectorXd& q, 
                         Eigen::VectorXd& p, 
                         Eigen::MatrixXd& phi, // State Transition Matrix
                         double ds);

    Eigen::Matrix3d compute_hessian(const Eigen::VectorXd& q) const;

    double compute_force_gradient(const Eigen::VectorXd& q) const;
    
    double compute_total_energy(const Eigen::VectorXd& q, const Eigen::VectorXd& p) const;

    void split_state(const Eigen::VectorXd& y,
                     Eigen::VectorXd& q, Eigen::VectorXd& p) const;
                     
    Eigen::VectorXd join_state(const Eigen::VectorXd& q,
                               const Eigen::VectorXd& p) const;

    double precision_;
    double mu_;
    double j2_;
    double r_eq_;
};

} // namespace astdyn::propagation

#endif // ASTDYN_AAS_INTEGRATOR_HPP
