#include "astdyn/propagation/OrbFitIntegrator.hpp"
#include <cmath>
#include <algorithm>
#include "astdyn/utils/Atomics.hpp"

namespace astdyn::propagation {

Eigen::VectorXd OrbFitDPIntegrator::integrate(const DerivativeFunction& f,
                                              const Eigen::VectorXd& y0,
                                              double t0,
                                              double tf) {
    std::vector<double> t_out;
    std::vector<Eigen::VectorXd> y_out;
    integrate_steps(f, y0, t0, tf, t_out, y_out);
    return y_out.back();
}

void OrbFitDPIntegrator::integrate_steps(const DerivativeFunction& f,
                                         const Eigen::VectorXd& y0,
                                         double t0,
                                         double tf,
                                         std::vector<double>& t_out,
                                         std::vector<Eigen::VectorXd>& y_out) {
    stats_.reset();
    
    // Dormand-Prince tableau
    static const double a21=1.0/5,
        a31=3.0/40, a32=9.0/40,
        a41=44.0/45, a42=-56.0/15, a43=32.0/9,
        a51=19372.0/6561, a52=-25360.0/2187, a53=64448.0/6561, a54=-212.0/729,
        a61=9017.0/3168, a62=-355.0/33, a63=46732.0/5247,
        a64=49.0/176, a65=-5103.0/18656;

    static const double b1=35.0/384, b3=500.0/1113, b4=125.0/192,
                        b5=-2187.0/6784, b6=11.0/84;
    static const double e1=71.0/57600, e3=-71.0/16695, e4=71.0/1920,
                        e5=-17253.0/339200, e6=22.0/525, e7=-1.0/40;

    Eigen::VectorXd y = y0;
    size_t dim = y0.size();
    double t = t0;
    double h = (tf > t) ? h_init_ : -h_init_;
    double sign = (tf > t) ? 1.0 : -1.0;

    t_out.push_back(t);
    y_out.push_back(y);

    int maxstep = 1000000;
    for(int step=0; step<maxstep; step++){
        if(sign*(tf - t) <= 0.0) break;
        if(std::abs(h) > std::abs(tf - t)) {
            h = tf - t;
        }

        Eigen::VectorXd k1 = f(t, y);
        Eigen::VectorXd k2 = f(t + h*0.2, y + h*a21*k1);
        Eigen::VectorXd k3 = f(t + h*0.3, y + h*(a31*k1 + a32*k2));
        Eigen::VectorXd k4 = f(t + h*0.8, y + h*(a41*k1 + a42*k2 + a43*k3));
        Eigen::VectorXd k5 = f(t + h,     y + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4));
        Eigen::VectorXd k6 = f(t + h,     y + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5));

        Eigen::VectorXd y5 = y + h*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6);
        Eigen::VectorXd k7 = f(t + h, y5);

        double err = 0.0;
        for(size_t i=0; i<dim; i++){
            double ei = h*(e1*k1[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*k7[i]);
            double sc = 1e-12 + tol_*std::max(std::abs(y[i]), std::abs(y5[i]));
            err += (ei/sc)*(ei/sc);
        }
        err = std::sqrt(err/dim);
        
        stats_.num_function_evals += 7;

        if(err <= 1.0){
            t += h;
            y  = y5;
            t_out.push_back(t);
            y_out.push_back(y);
            stats_.num_steps++;
            if(stats_.min_step_size == 0.0) {
                stats_.min_step_size = std::abs(h);
            } else {
                astdyn::utils::atomic_min(stats_.min_step_size, std::abs(h));
            }
            astdyn::utils::atomic_max(stats_.max_step_size, std::abs(h));
        } else {
            stats_.num_rejected_steps++;
        }

        double factor = 0.9 * std::pow(1.0/(err + 1e-10), 0.2);
        factor = std::max(0.2, std::min(factor, 10.0));
        h *= factor;
        h  = std::copysign(std::min(std::abs(h), h_max_), sign);
        h  = std::copysign(std::max(std::abs(h), h_min_), sign);
    }
    stats_.final_time = t;
}

std::vector<Eigen::VectorXd> OrbFitDPIntegrator::integrate_at(const DerivativeFunction& f,
                                                              const Eigen::VectorXd& y0,
                                                              double t0,
                                                              const std::vector<double>& t_targets) {
    std::vector<Eigen::VectorXd> results;
    results.reserve(t_targets.size());
    
    if (t_targets.empty()) return results;
    
    Eigen::VectorXd y = y0;
    double t = t0;
    
    for (double target_time : t_targets) {
        if (std::abs(target_time - t) < 1e-10) {
            results.push_back(y);
            continue;
        }
        y = integrate(f, y, t, target_time);
        results.push_back(y);
        t = target_time;
    }
    
    return results;
}

Eigen::VectorXd OrbFitRK4Integrator::integrate(const DerivativeFunction& f,
                                               const Eigen::VectorXd& y0,
                                               double t0,
                                               double tf) {
    std::vector<double> t_out;
    std::vector<Eigen::VectorXd> y_out;
    integrate_steps(f, y0, t0, tf, t_out, y_out);
    return y_out.back();
}

void OrbFitRK4Integrator::integrate_steps(const DerivativeFunction& f,
                                          const Eigen::VectorXd& y0,
                                          double t0,
                                          double tf,
                                          std::vector<double>& t_out,
                                          std::vector<Eigen::VectorXd>& y_out) {
    stats_.reset();
    
    Eigen::VectorXd y = y0;
    double t = t0;
    double h = (tf > t) ? std::abs(h_) : -std::abs(h_);
    double sign = (tf > t) ? 1.0 : -1.0;

    t_out.push_back(t);
    y_out.push_back(y);

    int maxstep = 1000000;
    for(int step=0; step<maxstep; step++){
        if(sign*(tf - t) <= 0.0) break;
        if(std::abs(h) > std::abs(tf - t)) {
            h = tf - t;
        }

        Eigen::VectorXd k1 = f(t, y);
        Eigen::VectorXd k2 = f(t + h/2.0, y + h*k1/2.0);
        Eigen::VectorXd k3 = f(t + h/2.0, y + h*k2/2.0);
        Eigen::VectorXd k4 = f(t + h,     y + h*k3);

        y = y + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        t += h;
        
        stats_.num_function_evals += 4;
        stats_.num_steps++;
        if(stats_.min_step_size == 0.0) {
            stats_.min_step_size = std::abs(h);
        } else {
            astdyn::utils::atomic_min(stats_.min_step_size, std::abs(h));
        }
        astdyn::utils::atomic_max(stats_.max_step_size, std::abs(h));

        t_out.push_back(t);
        y_out.push_back(y);
    }
    stats_.final_time = t;
}

std::vector<Eigen::VectorXd> OrbFitRK4Integrator::integrate_at(const DerivativeFunction& f,
                                                               const Eigen::VectorXd& y0,
                                                               double t0,
                                                               const std::vector<double>& t_targets) {
    std::vector<Eigen::VectorXd> results;
    results.reserve(t_targets.size());
    
    if (t_targets.empty()) return results;
    
    Eigen::VectorXd y = y0;
    double t = t0;
    
    for (double target_time : t_targets) {
        if (std::abs(target_time - t) < 1e-10) {
            results.push_back(y);
            continue;
        }
        y = integrate(f, y, t, target_time);
        results.push_back(y);
        t = target_time;
    }
    
    return results;
}


} // namespace astdyn::propagation
