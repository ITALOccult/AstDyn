/**
 * @file GRKNIntegrator.hpp
 * @brief GRKN — Generalized Runge-Kutta-Nyström integrator, embedded pair 6(4)
 *
 * Implements the RKN6(4)6FM formula of Dormand, El-Mikkawy & Prince (1987),
 * Table 5, IMA J. Numer. Anal. 7, 235–250.
 *
 * Key properties:
 *   - Order 6 for position (y),  order 4 for velocity (y')
 *   - 6 stages, FSAL (First Same As Last): 5 new NFE per accepted step
 *   - Designed for  y'' = f(x, y),  f independent of y'
 *   - Embedded 4th-order pair for adaptive step control
 */

#ifndef ASTDYN_GRKN_INTEGRATOR_HPP
#define ASTDYN_GRKN_INTEGRATOR_HPP

#include "astdyn/propagation/Integrator.hpp"
#include <Eigen/Dense>
#include <array>
#include <vector>

namespace astdyn::propagation {

static constexpr int GRKN_S = 6;   // number of stages, RKN6(4)6FM

class GRKNIntegrator : public Integrator {
public:
    /**
     * @param tolerance    Mixed abs/rel tolerance (default 1e-12)
     * @param initial_step Initial step size [days]
     * @param min_step     Minimum step [days]
     * @param max_step     Maximum step [days]
     * @param mu           GM central body [AU³/day²]  (GMS default)
     * @param j2           J2 coefficient  (0 = off)
     * @param r_eq         Equatorial radius for J2 [AU]
     */
    explicit GRKNIntegrator(
        double tolerance    = 1e-12,
        double initial_step = 1.0,
        double min_step     = 1e-8,
        double max_step     = 100.0,
        double mu           = 0.0002959122082855911,
        double j2           = 0.0,
        double r_eq         = 0.004650467261);

    Eigen::VectorXd integrate(
        const DerivativeFunction& f,
        const Eigen::VectorXd&    y0,
        double t0, double tf) override;

    void integrate_steps(
        const DerivativeFunction&     f,
        const Eigen::VectorXd&        y0,
        double t0, double tf,
        std::vector<double>&          t_out,
        std::vector<Eigen::VectorXd>& y_out) override;

    std::vector<Eigen::VectorXd> integrate_at(
        const DerivativeFunction&   f,
        const Eigen::VectorXd&      y0,
        double t0,
        const std::vector<double>&  t_targets) override;

    void set_tolerance(double tol)  { rtol_ = atol_ = tol; }
    void set_min_step(double h)     { h_min_ = h; }
    void set_max_step(double h)     { h_max_ = h; }
    void set_central_body(double mu, double j2 = 0.0, double r_eq = 4.6505e-3) {
        mu_ = mu; j2_ = j2; r_eq_ = r_eq;
    }

private:
    double rtol_, atol_;
    double h_initial_, h_min_, h_max_;
    double mu_, j2_, r_eq_;
    double err_prev_ = 1.0;

    // ── Butcher–Nyström tableau: RKN6(4)6FM  (DEP 1987, Table 5) ────────
    static constexpr std::array<double, GRKN_S> c_ = {{
        0.0,
        1.0  / 10.0,
        3.0  / 10.0,
        7.0  / 10.0,
        17.0 / 25.0,
        1.0
    }};

    static const double a_[GRKN_S][GRKN_S];

    static constexpr std::array<double, GRKN_S> bhat_ = {{
         151.0 /   2142.0,
           5.0 /    116.0,
         385.0 /   1368.0,
          55.0 /    168.0,
       -6250.0 /  28101.0,
           0.0
    }};

    static constexpr std::array<double, GRKN_S> bhatd_ = {{
          151.0 /   2142.0,
           25.0 /    522.0,
          275.0 /    684.0,
          275.0 /    252.0,
       -78125.0 / 112404.0,
            1.0 /     12.0
    }};

    static constexpr std::array<double, GRKN_S> b_ = {{
        1349.0 / 157500.0,
        7873.0 /  50000.0,
      192199.0 / 900000.0,
      521683.0 / 2100000.0,
         -16.0 /    125.0,
           0.0
    }};

    static constexpr std::array<double, GRKN_S> bd_ = {{
        1349.0 / 157500.0,
        7873.0 /  45000.0,
       27457.0 /  90000.0,
      521683.0 / 630000.0,
          -2.0 /      5.0,
           1.0 /     12.0
    }};

    using Vec3   = Eigen::Vector3d;
    using Mat3   = Eigen::Matrix3d;
    using Mat66  = Eigen::Matrix<double, 6, 6>;
    using Stages = std::array<Vec3, GRKN_S>;

    void rkn_step(
        const DerivativeFunction& f,
        double t, const Vec3& r, const Vec3& v, double h,
        Vec3& r6, Vec3& v6,
        Vec3& r4, Vec3& v4,
        Stages& K) const;

    void propagate_stm(
        const Vec3& r, const Vec3& v,
        double h, const Stages& K,
        Mat66& phi) const;

    Mat3 jacobian(const Vec3& r) const;

    Vec3 stage_pos(const Vec3& r, const Vec3& v, double h,
                   int i, const Stages& K) const;

    double error_norm(const Vec3& r_old,
                      const Vec3& r6, const Vec3& r4) const;

    double pi_controller(double err, double h) const;

    static void split(const Eigen::VectorXd& y, Vec3& r, Vec3& v);
    static Eigen::VectorXd join(const Vec3& r, const Vec3& v,
                                const Mat66* phi = nullptr);

    void adaptive_loop(
        const DerivativeFunction&     f,
        Vec3& r, Vec3& v, Mat66& phi, bool has_stm,
        double& t, double t_end,
        const std::vector<double>&    checkpoints,
        std::vector<double>&          t_out,
        std::vector<Eigen::VectorXd>& y_out);
};

} // namespace astdyn::propagation

#endif // ASTDYN_GRKN_INTEGRATOR_HPP
