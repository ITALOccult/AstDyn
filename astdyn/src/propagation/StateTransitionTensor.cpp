/**
 * @file StateTransitionTensor.cpp
 * @brief Joint (state, Phi, Psi) propagation, AAS adaptive Drift-Kick-Drift.
 *
 * Step control replicates AASIntegrator::estimate_step_size (force-gradient
 * metric, harmonic mean of the endpoint step lengths, min-period ceiling,
 * clamp).  The STM/STT sub-updates are Eqs. (B7)-(B11); the composition is
 * the seven-stage Yoshida-4 palindrome (Eq. B12).
 */

#include "astdyn/propagation/StateTransitionTensor.hpp"
#include <iostream>
#include <cstdlib>
#include "astdyn/core/Constants.hpp"
#include <algorithm>
#include <cmath>

namespace astdyn::propagation {

using astdyn::Vector3d;
using astdyn::Matrix3d;
using astdyn::Matrix6d;
using astdyn::math::Tensor3;
using astdyn::math::Tensor6;

StateTransitionTensor::StateTransitionTensor(PotentialModel model, double precision)
    : model_(std::move(model)), precision_(precision) {
    const double k = std::pow(2.0, 1.0 / 3.0);
    w1_ = 1.0 / (2.0 - k);
    w0_ = 1.0 - 2.0 * w1_;
    d1_ = w1_;  d2_ = w0_;
    c1_ = w1_ / 2.0;  c2_ = (w1_ + w0_) / 2.0;
}

// --- AAS force-gradient metric: max local |2 mu / r^3| (+ J2), all attractors -
double StateTransitionTensor::force_gradient(const Vector3d& r,
                                             const PotentialModel& m) const {
    const double rc = r.norm();
    double g = std::abs(2.0 * m.central_gm / (rc * rc * rc));
    if (m.j2 != 0.0) {
        const double r2 = rc * rc;
        g += std::abs(4.0 * m.central_gm * m.j2 * m.r_eq * m.r_eq / (r2 * r2 * rc));
    }
    for (size_t l = 0; l < m.perturber_gm.size(); ++l) {
        const double u = (r - m.perturber_pos[l]).norm();
        if (u > 1e-9) g = std::max(g, std::abs(2.0 * m.perturber_gm[l] / (u * u * u)));
    }
    return g;
}

double StateTransitionTensor::estimate_step(const Vector3d& q, const Vector3d& p,
                                            double target_dt, const PotentialModel& m) const {
    double min_r = q.norm();
    for (const auto& s : m.perturber_pos) min_r = std::min(min_r, (q - s).norm());
    if (!std::isfinite(min_r)) return 1e-8;

    const double sp = precision_ * 50.0;
    const double g_n = sp / std::sqrt(std::max(force_gradient(q, m), 1e-20));
    const Vector3d q_next = q + p * g_n;
    const double g_p = sp / std::sqrt(std::max(force_gradient(q_next, m), 1e-20));
    double g_avg = 2.0 * g_n * g_p / (g_n + g_p);

    const double T_min = astdyn::constants::TWO_PI
                       * std::sqrt(min_r * min_r * min_r / m.central_gm);
    g_avg = std::min(g_avg, T_min / 20.0);
    g_avg = std::clamp(g_avg, 1e-8, 100.0);
    return std::min(g_avg, target_dt);
}

// --- Kick K(h): p += h a; momentum-row STM/STT updates (Eqs. B7-B8) -----------
void StateTransitionTensor::kick(const Vector3d& q, Vector3d& p,
                                 Matrix6d& phi, Tensor6& psi,
                                 double h, const PotentialModel& m) const {
    const Matrix3d G = -potential_hessian(q, m);          // acc. gradient = -U_ij
    Tensor3 T3 = potential_third_derivative(q, m);
    T3 *= -1.0;                                           // = -U_ijk
    p += h * acceleration(q, m);

    const Eigen::Matrix<double, 3, 6> Phi_pos = phi.topRows<3>();
    phi.bottomRows<3>() += h * (G * Phi_pos);                        // (B7)
    for (int i = 0; i < 3; ++i) {                                   // (B8)
        Matrix6d gpsi = Matrix6d::Zero();
        for (int j = 0; j < 3; ++j) gpsi += G(i, j) * psi.slice(j);
        const Matrix6d src = T3.contract_bilinear<6>(i, Phi_pos);
        psi.slice(3 + i) += h * (gpsi + src);
    }
}

// --- Drift D(h): q += h p; position-row STM/STT updates (Eqs. B10-B11) --------
void StateTransitionTensor::drift(Vector3d& q, const Vector3d& p,
                                  Matrix6d& phi, Tensor6& psi, double h) const {
    q += h * p;
    phi.topRows<3>() += h * phi.bottomRows<3>();                    // (B10)
    for (int i = 0; i < 3; ++i)
        psi.slice(i) += h * psi.slice(3 + i);                      // (B11)
}

// --- one seven-stage Yoshida-4 macro-step (Eq. B12) --------------------------
void StateTransitionTensor::step(Vector3d& q, Vector3d& p,
                                 Matrix6d& phi, Tensor6& psi,
                                 double dt, const PotentialModel& m) const {
    drift(q, p, phi, psi, c1_ * dt);
    kick (q, p, phi, psi, d1_ * dt, m);
    drift(q, p, phi, psi, c2_ * dt);
    kick (q, p, phi, psi, d2_ * dt, m);
    drift(q, p, phi, psi, c2_ * dt);
    kick (q, p, phi, psi, d1_ * dt, m);
    drift(q, p, phi, psi, c1_ * dt);
}

StateTransitionTensor::State
StateTransitionTensor::to_state(const Vector3d& q, const Vector3d& p,
                               time::EpochTDB t,
                               const astdyn::physics::GravitationalParameter& gm) const {
    using astdyn::physics::Distance;
    using astdyn::physics::Velocity;
    const auto pos = astdyn::math::Vector3<Frame, Distance>::from_si(
        Distance::from_au(q[0]).to_m(), Distance::from_au(q[1]).to_m(),
        Distance::from_au(q[2]).to_m());
    const auto vel = astdyn::math::Vector3<Frame, Velocity>::from_si(
        Velocity::from_au_d(p[0]).to_ms(), Velocity::from_au_d(p[1]).to_ms(),
        Velocity::from_au_d(p[2]).to_ms());
    return State(t, pos, vel, gm);
}

StateTransitionTensor::Result
StateTransitionTensor::propagate(const State& x0, time::EpochTDB t1,
                                 const EphemerisRefresh& refresh) const {
    PotentialModel m = model_;                                  // per-call, refreshable
    const Eigen::Matrix<double, 6, 1> y0 = x0.to_eigen_au_aud();
    Vector3d q = y0.head<3>();
    Vector3d p = y0.tail<3>();

    Matrix6d phi = Matrix6d::Identity();
    Tensor6 psi;                                            // zero-initialised

    const double span = (t1 - x0.epoch).to_days();          // forward: span >= 0
    const double total = std::abs(span);
    const double dir = (span >= 0.0) ? 1.0 : -1.0;
    double elapsed = 0.0;
    // Il numero di passi e la loro ampiezza: il dato che non ho mai guardato
    // mentre tiravo a indovinare sulla causa dell'errore.
    long n_steps = 0;
    double g_min = 1e300, g_max = 0.0;

    while (elapsed < total - 1e-12) {
        const auto t_now = x0.epoch + time::TimeDuration::from_days(dir * elapsed);

        // Il passo si stima con i perturbatori dove sono adesso.
        if (refresh) refresh(t_now, m);
        const double g = estimate_step(q, p, total - elapsed, m);

        // Ma il passo li tratta come FERMI per tutta la sua durata, e in dieci
        // giorni Giove fa undici milioni di km. Valutarli all'inizio del passo
        // e' un errore del primo ordine, e per giunta sistematico: sbaglia
        // sempre nello stesso verso e si accumula invece di mediarsi.
        // Valutarli a META' PASSO lo rende del secondo ordine, al prezzo di una
        // seconda chiamata.
        //
        // Resta un'approssimazione: un passo esatto vorrebbe i perturbatori a
        // ogni stadio dell'integratore, non uno solo per passo.
        if (refresh) refresh(t_now + time::TimeDuration::from_days(dir * 0.5 * g), m);

        step(q, p, phi, psi, dir * g, m);
        elapsed += g;
        ++n_steps;
        g_min = std::min(g_min, g);
        g_max = std::max(g_max, g);
    }
    if (std::getenv("ASTDYN_STT_STEPS")) {
        std::cerr << "[stt] " << n_steps << " passi su " << total << " giorni"
                  << "   g in [" << g_min << ", " << g_max << "] d"
                  << "   precisione " << precision_ << "\n";
    }
    return Result{to_state(q, p, t1, x0.gm), phi, psi};
}

StateTransitionTensor::Result
StateTransitionTensor::propagate_fixed(const State& x0, time::EpochTDB t1,
                                       int n_steps) const {
    PotentialModel m = model_;
    const Eigen::Matrix<double, 6, 1> y0 = x0.to_eigen_au_aud();
    Vector3d q = y0.head<3>();
    Vector3d p = y0.tail<3>();

    Matrix6d phi = Matrix6d::Identity();
    Tensor6 psi;

    const double dt = (t1 - x0.epoch).to_days() / static_cast<double>(n_steps);
    for (int i = 0; i < n_steps; ++i) step(q, p, phi, psi, dt, m);
    return Result{to_state(q, p, t1, x0.gm), phi, psi};
}

}  // namespace astdyn::propagation
