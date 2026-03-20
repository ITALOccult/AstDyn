#include "astdyn/propagation/GRKNIntegrator.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace astdyn::propagation {

const double GRKNIntegrator::a_[GRKN_S][GRKN_S] = {
    {0, 0, 0, 0, 0, 0},                                        // Row 0
    {1.0 / 200.0, 0, 0, 0, 0, 0},                              // Row 1
    {-1.0 / 2200.0, 1.0 / 22.0, 0, 0, 0, 0},                   // Row 2
    {637.0 / 6600.0, -7.0 / 110.0, 7.0 / 33.0, 0, 0, 0},       // Row 3
    {225437.0 / 1968750.0, -30073.0 / 281250.0, 65569.0 / 281250.0, -9367.0 / 984375.0, 0, 0}, // Row 4
    {151.0 / 2142.0, 5.0 / 116.0, 385.0 / 1368.0, 55.0 / 168.0, -6250.0 / 28101.0, 0}         // Row 5
};

GRKNIntegrator::GRKNIntegrator(double tolerance, double initial_step, double min_step, double max_step, double mu, double j2, double r_eq)
    : rtol_(tolerance), atol_(tolerance), h_initial_(initial_step), h_min_(min_step), h_max_(max_step), mu_(mu), j2_(j2), r_eq_(r_eq) {}

void GRKNIntegrator::split(const Eigen::VectorXd& y, Vec3& r, Vec3& v) {
    r = y.head<3>();
    v = y.segment<3>(3);
}

Eigen::VectorXd GRKNIntegrator::join(const Vec3& r, const Vec3& v, const Mat66* phi) {
    if (phi) {
        Eigen::VectorXd y(42);
        y.head<3>() = r;
        y.segment<3>(3) = v;
        Eigen::Map<Eigen::MatrixXd>(y.data() + 6, 6, 6) = *phi;
        return y;
    } else {
        Eigen::VectorXd y(6);
        y.head<3>() = r;
        y.segment<3>(3) = v;
        return y;
    }
}

double GRKNIntegrator::error_norm(const Vec3& r_old, const Vec3& r6, const Vec3& r4) const {
    double sum = 0;
    for (int i = 0; i < 3; ++i) {
        double sc = atol_ + rtol_ * std::max(std::abs(r_old[i]), std::abs(r6[i]));
        double e = (r6[i] - r4[i]) / sc;
        sum += e * e;
    }
    return std::sqrt(sum / 3.0);
}

double GRKNIntegrator::pi_controller(double err, double h) const {
    const double safety = 0.9;
    const double p = 4.0; // embedded order
    const double alpha = 0.7 / p;
    const double beta = 0.4 / p;
    
    double fac = safety * std::pow(std::max(err, 1e-18), -alpha) * std::pow(std::max(err_prev_, 1e-18), beta);
    fac = std::clamp(fac, 0.1, 5.0);
    return h * fac;
}

GRKNIntegrator::Vec3 GRKNIntegrator::stage_pos(const Vec3& r, const Vec3& v, double h, int i, const Stages& K) const {
    Vec3 Ri = r + c_[i] * h * v;
    double h2 = h * h;
    for (int j = 0; j < i; ++j) {
        Ri += h2 * a_[i][j] * K[j];
    }
    return Ri;
}

void GRKNIntegrator::rkn_step(const DerivativeFunction& f, double t, const Vec3& r, const Vec3& v, double h, Vec3& r6, Vec3& v6, Vec3& r4, Vec3& v4, Stages& K) const {
    // Stage evaluation
    for (int i = 0; i < GRKN_S; ++i) {
        Vec3 Ri = stage_pos(r, v, h, i, K);
        // Specialised force evaluation: velocity is zero
        Eigen::VectorXd y_dummy = Eigen::VectorXd::Zero(6);
        y_dummy.head<3>() = Ri;
        K[i] = f(t + c_[i] * h, y_dummy).tail<3>();
    }
    
    // Position/Velocity updates
    double h2 = h * h;
    r6 = r + h * v;
    v6 = v;
    r4 = r + h * v;
    v4 = v;
    
    for (int i = 0; i < GRKN_S; ++i) {
        r6 += h2 * bhat_[i] * K[i];
        v6 += h * bhatd_[i] * K[i];
        r4 += h2 * b_[i] * K[i];
        v4 += h * bd_[i] * K[i];
    }
}

GRKNIntegrator::Mat3 GRKNIntegrator::jacobian(const Vec3& q) const {
    const double r = q.norm();
    const double r2 = r * r;
    const double r3 = r2 * r;
    const double r5 = r3 * r2;
    
    Mat3 I = Mat3::Identity();
    Mat3 H = (-mu_ / r3) * I + (3.0 * mu_ / r5) * (q * q.transpose());
    
    if (j2_ != 0.0) {
        const double req2 = r_eq_ * r_eq_;
        const double z = q[2];
        const double z2 = z * z;
        const double r7 = r5 * r2;
        const double r9 = r7 * r2;
        const double C = 0.5 * mu_ * j2_ * req2;
        
        Mat3 H_j2;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double delta_ij = (i == j) ? 1.0 : 0.0;
                double delta_iz = (i == 2) ? 1.0 : 0.0;
                double delta_jz = (j == 2) ? 1.0 : 0.0;
                double term1 = (15.0 * z2 / r7 - 3.0 / r5) * delta_ij;
                double term2 = q[i] * (30.0 * z * delta_jz / r7 - 105.0 * z2 * q[j] / r9 + 15.0 * q[j] / r7);
                double term3 = -6.0 * delta_iz * (delta_jz / r5 - 5.0 * z * q[j] / r7);
                H_j2(i, j) = C * (term1 + term2 + term3);
            }
        }
        H += H_j2;
    }
    return H;
}

void GRKNIntegrator::propagate_stm(const Vec3& r, const Vec3& v, double h, const Stages& K, Mat66& phi) const {
    // δr'' = J(r) δr
    // We update phi using the Nyström formulation on the variational equations
    std::array<Mat3, GRKN_S> J_i;
    for (int i = 0; i < GRKN_S; ++i) {
        J_i[i] = jacobian(stage_pos(r, v, h, i, K));
    }
    
    // Stage variations
    std::array<Eigen::Matrix<double, 3, 6>, GRKN_S> dR;
    auto dR0 = phi.block<3, 6>(0, 0); // initial pos variation
    auto dV0 = phi.block<3, 6>(3, 0); // initial vel variation
    
    double h2 = h * h;
    for (int i = 0; i < GRKN_S; ++i) {
        dR[i] = dR0 + c_[i] * h * dV0;
        for (int j = 0; j < i; ++j) {
            dR[i] += h2 * a_[i][j] * J_i[j] * dR[j];
        }
    }
    
    // Final update
    phi.block<3, 6>(0, 0) += h * dV0; // base linear part
    for (int i = 0; i < GRKN_S; ++i) {
        phi.block<3, 6>(0, 0) += h2 * bhat_[i] * J_i[i] * dR[i];
        phi.block<3, 6>(3, 0) += h * bhatd_[i] * J_i[i] * dR[i];
    }
}

void GRKNIntegrator::adaptive_loop(const DerivativeFunction& f, Vec3& r, Vec3& v, Mat66& phi, bool has_stm,
                                  double& t, double t_end, const std::vector<double>& checkpoints,
                                  std::vector<double>& t_out, std::vector<Eigen::VectorXd>& y_out) {
    double h = h_initial_ * ((t_end > t) ? 1.0 : -1.0);
    int cp_idx = 0;
    Stages K;
    err_prev_ = 1.0;
    bool first_step = true;

    while (std::abs(t_end - t) > 1e-14) {
        if (std::abs(h) > std::abs(t_end - t)) h = t_end - t;
        
        Vec3 r6, v6, r4, v4;
        rkn_step(f, t, r, v, h, r6, v6, r4, v4, K);
        
        if (!r6.allFinite() || !v6.allFinite()) {
            stats_.num_rejected_steps++;
            h *= 0.5;
            if (std::abs(h) < h_min_) break;
            continue;
        }

        double err = error_norm(r, r6, r4);
        if (err <= 1.0 || std::abs(h) <= h_min_) {
            // Accepted
            if (has_stm) propagate_stm(r, v, h, K, phi);
            
            r = r6; v = v6; t += h;
            stats_.num_steps++;
            err_prev_ = err;
            
            while (cp_idx < static_cast<int>(checkpoints.size()) && 
                   ((h > 0 && checkpoints[cp_idx] <= t) || (h < 0 && checkpoints[cp_idx] >= t))) {
                t_out.push_back(checkpoints[cp_idx]);
                y_out.push_back(join(r, v, has_stm ? &phi : nullptr));
                cp_idx++;
            }
            
            h = pi_controller(err, h);
            h = std::clamp(std::abs(h), h_min_, h_max_) * ((t_end > t) ? 1.0 : -1.0);
        } else {
            // Rejected
            stats_.num_rejected_steps++;
            h *= 0.5;
            if (std::abs(h) < h_min_) h = h_min_ * ((t_end > t) ? 1.0 : -1.0);
        }
        
        if (stats_.num_steps > 1000000) break;
    }
}

Eigen::VectorXd GRKNIntegrator::integrate(const DerivativeFunction& f, const Eigen::VectorXd& y0, double t0, double tf) {
    stats_.reset();
    Vec3 r, v; split(y0, r, v);
    Mat66 phi; bool has_stm = (y0.size() == 42);
    if (has_stm) phi = Eigen::Map<const Mat66>(y0.data() + 6);
    
    double t = t0;
    std::vector<double> t_out; std::vector<Eigen::VectorXd> y_out;
    adaptive_loop(f, r, v, phi, has_stm, t, tf, {}, t_out, y_out);
    
    return join(r, v, has_stm ? &phi : nullptr);
}

void GRKNIntegrator::integrate_steps(const DerivativeFunction& f, const Eigen::VectorXd& y0, double t0, double tf, std::vector<double>& t_out, std::vector<Eigen::VectorXd>& y_out) {
    stats_.reset();
    Vec3 r, v; split(y0, r, v);
    Mat66 phi; bool has_stm = (y0.size() == 42);
    if (has_stm) phi = Eigen::Map<const Mat66>(y0.data() + 6);
    
    t_out.push_back(t0); y_out.push_back(y0);
    double t = t0;
    adaptive_loop(f, r, v, phi, has_stm, t, tf, {}, t_out, y_out);
}

std::vector<Eigen::VectorXd> GRKNIntegrator::integrate_at(const DerivativeFunction& f, const Eigen::VectorXd& y0, double t0, const std::vector<double>& t_targets) {
    stats_.reset();
    Vec3 r, v; split(y0, r, v);
    Mat66 phi; bool has_stm = (y0.size() == 42);
    if (has_stm) phi = Eigen::Map<const Mat66>(y0.data() + 6);
    
    std::vector<double> t_out; std::vector<Eigen::VectorXd> res;
    double t = t0;
    adaptive_loop(f, r, v, phi, has_stm, t, t_targets.back(), t_targets, t_out, res);
    return res;
}

} // namespace astdyn::propagation
