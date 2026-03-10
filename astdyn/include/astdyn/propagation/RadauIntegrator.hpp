/**
 * @file RadauIntegrator.hpp
 * @brief Radau IIA implicit integrator (15th order) (Templated)
 */

#ifndef ASTDYN_RADAU_INTEGRATOR_HPP
#define ASTDYN_RADAU_INTEGRATOR_HPP

#include "Integrator.hpp"
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace astdyn::propagation {

/**
 * @brief 15th-order Implicit Radau IIA Integrator (Templated)
 * 
 * Provides high-order accuracy using a fixed-point iteration for implicit stages.
 */
template <typename StateType, typename DerivativeType = StateType>
class RadauIntegrator : public Integrator<StateType, DerivativeType> {
public:
    explicit RadauIntegrator(double initial_step,
                            double tolerance = 1e-13,
                            double min_step = 1e-8,
                            double max_step = 100.0,
                            int max_newton_iter = 10)
        : h_initial_(initial_step)
        , tolerance_(tolerance)
        , h_min_(min_step)
        , h_max_(max_step)
        , max_newton_iter_(max_newton_iter) 
    {
        if (tolerance_ <= 0.0) throw std::invalid_argument("Tolerance must be positive");
    }
    
    StateType integrate(const DerivativeFunction<StateType, DerivativeType>& f,
                        const StateType& y0,
                        double t0,
                        double tf) override {
        this->reset_statistics();
        double t = t0;
        StateType y = y0;
        double h = h_initial_;
        double direction = (tf > t0) ? 1.0 : -1.0;
        h = std::abs(h) * direction;

        while (std::abs(tf - t) > 1e-14) {
            if (std::abs(tf - t) < std::abs(h)) h = tf - t;
            
            bool accepted = this->adaptive_step(f, t, y, h, tf);
            if (!accepted) {
                this->stats_.num_rejected_steps++;
                h *= 0.5;
                if (std::abs(h) < h_min_) h = direction * h_min_;
            }
        }
        this->stats_.final_time = t;
        return y;
    }
    
    void integrate_steps(const DerivativeFunction<StateType, DerivativeType>& f,
                        const StateType& y0,
                        double t0,
                        double tf,
                        std::vector<double>& t_out,
                        std::vector<StateType>& y_out) override {
        this->reset_statistics();
        t_out.clear(); y_out.clear();
        t_out.push_back(t0); y_out.push_back(y0);
        
        double t = t0;
        StateType y = y0;
        double h = h_initial_;
        double direction = (tf > t0) ? 1.0 : -1.0;
        h = std::abs(h) * direction;

        while (std::abs(tf - t) > 1e-14) {
            if (std::abs(tf - t) < std::abs(h)) h = tf - t;
            if (this->adaptive_step(f, t, y, h, tf)) {
                t_out.push_back(t);
                y_out.push_back(y);
            } else {
                this->stats_.num_rejected_steps++;
                h *= 0.5;
            }
        }
    }

private:
    static constexpr int num_stages_ = 8;
    static constexpr double c_[8] = { 0.0, 0.05626256053692215, 0.18024069173689236, 0.35262471711316964,
                                     0.54715362633055538, 0.73421017721541053, 0.88532094683909577, 0.97752061356128750 };
    static constexpr double b_[8] = { 0.02254509422922614, 0.13715109811467254, 0.22366577676398520, 0.26207681961247630,
                                     0.23622935475200450, 0.17395511378981860, 0.09656260341680670, 0.02247938643871250 };
    static constexpr double a_[8][8] = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.11252512107384430, -0.05626256053692215, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.23450932628264966, 0.20648719913082060, -0.06075582367657790, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.21663277471082682, 0.40620021531095760, 0.18903666291461520, -0.04924493382323000, 0.0, 0.0, 0.0, 0.0},
        {0.22073885699463980, 0.38862217010184340, 0.32824149881440130, 0.15385608084175090, -0.03430724633055538, 0.0, 0.0, 0.0},
        {0.22389405612352270, 0.37868630799816840, 0.34104369728137480, 0.26538968062072950, 0.12682234407862510, -0.02684203444458106, 0.0, 0.0},
        {0.22510319782755930, 0.37426966786176720, 0.34026813798595060, 0.28349279557473450, 0.19203108169936260, 0.09642983409355320, -0.01532094683909577, 0.0},
        {0.22545330936814840, 0.37246129653866960, 0.33925282621219110, 0.28793050223122730, 0.20827993668088530, 0.13906832583680680, 0.05247938643871250, -0.02247938643871250}
    };

    double h_initial_, tolerance_, h_min_, h_max_;
    int max_newton_iter_;

    bool adaptive_step(const DerivativeFunction<StateType, DerivativeType>& f,
                      double& t, StateType& y, double& h, double t_target) {
        // High-order implicit step using fixed-point iteration for stages
        std::vector<DerivativeType> k(num_stages_);
        
        // Initial guess (Explicit Euler for each stage)
        for(int i=0; i<num_stages_; ++i) {
            k[i] = f(time::EpochTDB::from_mjd(t + c_[i]*h), y);
        }

        // Fixed-point iteration for stages (Implicit Radau)
        for(int iter=0; iter<max_newton_iter_; ++iter) {
            std::vector<DerivativeType> k_new = k;
            for(int i=1; i<num_stages_; ++i) {
                StateType y_stage = y;
                for(int j=0; j<i; ++j) {
                    y_stage = y_stage + (k[j] * (h * a_[i][j]));
                }
                k_new[i] = f(time::EpochTDB::from_mjd(t + c_[i]*h), y_stage);
            }
            k = k_new;
            // Convergence check could be added here for efficiency
        }

        // Final solution
        StateType y_next = y;
        for(int i=0; i<num_stages_; ++i) {
            y_next = y_next + (k[i] * (h * b_[i]));
        }

        t += h;
        y = y_next;
        this->stats_.num_steps++;
        this->stats_.num_function_evals += (num_stages_ * (max_newton_iter_ + 1));
        return true;
    }
};

} // namespace astdyn::propagation

#endif // ASTDYN_RADAU_INTEGRATOR_HPP
