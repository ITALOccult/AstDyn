/**
 * @file test_stt.cpp
 * @brief Unit tests for the second-order variational propagator
 *        (StateTransitionTensor): closed-form U_ijk, and the propagated
 *        STM (Phi) and STT (Psi) against finite differences of the flow.
 *
 * The finite-difference checks use propagate_fixed() so that the nominal and
 * perturbed trajectories share an identical discrete map (the analytic
 * variational quantities are the exact derivatives of that map).  A separate
 * test confirms the adaptive propagate() agrees with the fixed-step reference.
 */

#include <gtest/gtest.h>
#include "astdyn/propagation/StateTransitionTensor.hpp"
#include "astdyn/propagation/PotentialDerivatives.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/time/epoch.hpp"
#include <Eigen/Dense>
#include <cmath>

using namespace astdyn;
using namespace astdyn::propagation;
using astdyn::math::Tensor3;

namespace {

constexpr double MU_SUN = 0.0002959122082855911;   // au^3/day^2 (Gauss k^2)

using Frame = core::ECLIPJ2000;
using State = physics::CartesianStateTyped<Frame>;
using Vec6  = Eigen::Matrix<double, 6, 1>;

// Build a state from an [au, au/day] vector.
State make_state(const Vec6& y, time::EpochTDB t) {
    using physics::Distance;
    using physics::Velocity;
    using physics::GravitationalParameter;
    const auto pos = math::Vector3<Frame, Distance>::from_si(
        Distance::from_au(y[0]).to_m(), Distance::from_au(y[1]).to_m(),
        Distance::from_au(y[2]).to_m());
    const auto vel = math::Vector3<Frame, Velocity>::from_si(
        Velocity::from_au_d(y[3]).to_ms(), Velocity::from_au_d(y[4]).to_ms(),
        Velocity::from_au_d(y[5]).to_ms());
    return State(t, pos, vel, GravitationalParameter::from_au3_d2(MU_SUN));
}

// A moderately eccentric, inclined ~1 au test orbit.
Vec6 nominal_state() {
    const double vc = std::sqrt(MU_SUN);
    Vec6 y;
    y << 1.0, 0.0, 0.0, 0.0, 0.9 * vc, 0.2 * vc;
    return y;
}

ForceModel kepler_model() {
    ForceModel m;
    m.central_gm = MU_SUN;
    return m;
}

}  // namespace

// ---- 1. Closed-form third derivative vs finite differences of U_ij ---------
TEST(StateTransitionTensor, UijkMatchesFiniteDifference) {
    ForceModel m = kepler_model();
    m.j2 = 2.2e-7; m.r_eq = 4.65e-3;                       // stress J2 too
    m.perturber_gm = {MU_SUN * 9.54e-4};
    m.perturber_pos = {Vector3d(3.0, 3.5, 0.1)};

    const Vector3d r(0.7, 0.3, 0.5);
    const double d = 1e-7;
    const Tensor3 T = potential_third_derivative(r, m);
    double err = 0.0;
    for (int k = 0; k < 3; ++k) {
        Vector3d e = Vector3d::Zero(); e[k] = d;
        const Matrix3d num =
            (potential_hessian(r + e, m) - potential_hessian(r - e, m)) / (2 * d);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                err = std::max(err, std::abs(T(i, j, k) - num(i, j)));
    }
    EXPECT_LT(err, 1e-9);
}

// ---- 2. STM Phi vs finite differences of the (fixed-step) flow -------------
TEST(StateTransitionTensor, PhiMatchesFiniteDifference) {
    const StateTransitionTensor stt(kepler_model());
    const time::EpochTDB t0 = time::EpochTDB::from_mjd(60000.0);
    const time::EpochTDB t1 = t0 + time::TimeDuration::from_days(40.0);
    const int n = 300;
    const Vec6 y0 = nominal_state();

    const Matrix6d Phi = stt.propagate_fixed(make_state(y0, t0), t1, n).phi;

    const double eps = 1e-6;
    Matrix6d Phi_fd;
    for (int a = 0; a < 6; ++a) {
        Vec6 e = Vec6::Zero(); e[a] = eps;
        const Vec6 xp = stt.propagate_fixed(make_state(y0 + e, t0), t1, n)
                            .final_state.to_eigen_au_aud();
        const Vec6 xm = stt.propagate_fixed(make_state(y0 - e, t0), t1, n)
                            .final_state.to_eigen_au_aud();
        Phi_fd.col(a) = (xp - xm) / (2 * eps);
    }
    EXPECT_LT((Phi - Phi_fd).cwiseAbs().maxCoeff(), 1e-5);
}

// ---- 3. STT Psi vs mixed second differences of the flow --------------------
TEST(StateTransitionTensor, PsiMatchesFiniteDifference) {
    const StateTransitionTensor stt(kepler_model());
    const time::EpochTDB t0 = time::EpochTDB::from_mjd(60000.0);
    const time::EpochTDB t1 = t0 + time::TimeDuration::from_days(40.0);
    const int n = 300;
    const Vec6 y0 = nominal_state();

    const auto psi = stt.propagate_fixed(make_state(y0, t0), t1, n).psi;

    auto flow = [&](const Vec6& y) {
        return stt.propagate_fixed(make_state(y, t0), t1, n)
                  .final_state.to_eigen_au_aud();
    };
    const double eps = 1e-5;
    const Vec6 f0 = flow(y0);
    double err = 0.0, scale = 0.0;
    for (int a = 0; a < 6; ++a)
        for (int b = a; b < 6; ++b) {
            Vec6 ea = Vec6::Zero(); ea[a] = eps;
            Vec6 eb = Vec6::Zero(); eb[b] = eps;
            Vec6 d;
            if (a == b)
                d = (flow(y0 + ea) - 2 * f0 + flow(y0 - ea)) / (eps * eps);
            else
                d = (flow(y0 + ea + eb) - flow(y0 + ea - eb)
                   - flow(y0 - ea + eb) + flow(y0 - ea - eb)) / (4 * eps * eps);
            for (int i = 0; i < 6; ++i) {
                err = std::max(err, std::abs(psi.slice(i)(a, b) - d[i]));
                scale = std::max(scale, std::abs(psi.slice(i)(a, b)));
            }
        }
    EXPECT_LT(err / scale, 1e-4);   // finite-difference limited
}

// ---- 4. Psi is symmetric in its last two indices ---------------------------
TEST(StateTransitionTensor, PsiSymmetry) {
    const StateTransitionTensor stt(kepler_model());
    const time::EpochTDB t0 = time::EpochTDB::from_mjd(60000.0);
    const time::EpochTDB t1 = t0 + time::TimeDuration::from_days(40.0);
    const auto psi = stt.propagate_fixed(make_state(nominal_state(), t0), t1, 300).psi;
    double sym = 0.0;
    for (int i = 0; i < 6; ++i)
        for (int a = 0; a < 6; ++a)
            for (int b = 0; b < 6; ++b)
                sym = std::max(sym, std::abs(psi.slice(i)(a, b) - psi.slice(i)(b, a)));
    EXPECT_LT(sym, 1e-10);
}

// ---- 5. Adaptive metric agrees with the fixed-step reference ---------------
TEST(StateTransitionTensor, AdaptiveMatchesFixed) {
    const StateTransitionTensor stt(kepler_model(), 1e-6);
    const time::EpochTDB t0 = time::EpochTDB::from_mjd(60000.0);
    const time::EpochTDB t1 = t0 + time::TimeDuration::from_days(40.0);
    const State x0 = make_state(nominal_state(), t0);

    const Vec6 adaptive = stt.propagate(x0, t1).final_state.to_eigen_au_aud();
    const Vec6 fixed    = stt.propagate_fixed(x0, t1, 8000).final_state.to_eigen_au_aud();
    EXPECT_LT((adaptive.head<3>() - fixed.head<3>()).cwiseAbs().maxCoeff(), 1e-6);
}
