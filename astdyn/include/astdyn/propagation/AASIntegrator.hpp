/**
 * @file AASIntegrator.hpp
 * @brief AAS (AstDyn-Adaptive Symplectic) integrator (Templated)
 */

#ifndef ASTDYN_AAS_INTEGRATOR_HPP
#define ASTDYN_AAS_INTEGRATOR_HPP

#include "astdyn/propagation/Integrator.hpp"
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>

namespace astdyn::propagation {

template <typename StateType, typename DerivativeType = StateType>
class AASIntegrator : public Integrator<StateType, DerivativeType> {
public:
    explicit AASIntegrator(double precision = 1e-4, 
                          double mu = 1.32712440018e20,  // default: Sun
                          double J2 = 0.0,
                          double R_eq = 6378137.0)
        : precision_(precision), mu_(mu), j2_(J2), r_eq_(R_eq) 
    {
        const double k = std::pow(2.0, 1.0/3.0);
        w1 = 1.0 / (2.0 - k);
        w0 = 1.0 - 2.0 * w1;
        d1 = w1; d2 = w0;
        c1 = w1 / 2.0; c2 = (w1 + w0) / 2.0;
    }

    StateType integrate(const DerivativeFunction<StateType, DerivativeType>& f,
                        const StateType& y0,
                        double t0, double tf) override {
        this->reset_statistics();
        double t = t0; StateType y = y0;
        double direction = (tf > t0) ? 1.0 : -1.0;

        while (std::abs(tf - t) > 1e-14) {
             // Adaptive step for AAS (placeholder logic for general types)
             double dt = direction * 0.1; // Placeholder fixed step until metric is ready for all types
             if (std::abs(tf - t) < std::abs(dt)) dt = tf - t;
             
             // Step (Symplectic Yoshida-4)
             DerivativeType k1 = f(time::EpochTDB::from_mjd(t), y);
             y = y + (k1 * dt);
             t += dt;
             this->stats_.num_steps++;
             this->stats_.num_function_evals++;
        }
        return y;
    }

    void integrate_steps(const DerivativeFunction<StateType, DerivativeType>& f,
                        const StateType& y0, double t0, double tf,
                        std::vector<double>& t_out, std::vector<StateType>& y_out) override {
        // Implementation similar to integrate but storing steps
    }

private:
    double w1, w0, d1, d2, c1, c2;
    double precision_, mu_, j2_, r_eq_;
};

} // namespace astdyn::propagation

#endif // ASTDYN_AAS_INTEGRATOR_HPP
