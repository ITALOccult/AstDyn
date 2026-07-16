/**
 * @file PotentialDerivatives.cpp
 * @brief Closed-form U_ij and U_ijk (Keplerian + solar J2 + N-body).
 *
 * The J2 forms reproduce the AAS Hessian (Eq. A20) on one contraction; the
 * third-derivative closed forms are Eqs. (B1)-(B3) of the occultation paper,
 * cross-checked numerically against finite differences of U_ij.
 */

#include "astdyn/propagation/PotentialDerivatives.hpp"

namespace astdyn::propagation {

using astdyn::Vector3d;
using astdyn::Matrix3d;
using astdyn::math::Tensor3;

// ---- Keplerian U = -mu/r ------------------------------------------------
Matrix3d kepler_hessian(const Vector3d& u, double mu) {
    const double r = u.norm();
    const Vector3d rh = u / r;
    // U_ij = (mu/r^3)(delta_ij - 3 rh_i rh_j)
    return (mu / (r * r * r)) * (Matrix3d::Identity() - 3.0 * rh * rh.transpose());
}

Tensor3 kepler_third(const Vector3d& u, double mu) {
    const double r = u.norm();
    const Vector3d rh = u / r;
    const double c = 3.0 * mu / (r * r * r * r);
    Tensor3 T;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) {
                const double s1 = (i == j) * rh[k] + (i == k) * rh[j] + (j == k) * rh[i];
                T(i, j, k) = c * (5.0 * rh[i] * rh[j] * rh[k] - s1);
            }
    return T;
}

// ---- Solar J2 (pole = +z); C = 1/2 mu J2 R_eq^2 -------------------------
static double j2_C(double mu, double j2, double r_eq) {
    return 0.5 * mu * j2 * r_eq * r_eq;
}

Matrix3d j2_hessian(const Vector3d& r, double mu, double j2, double r_eq) {
    const double C = j2_C(mu, j2, r_eq);
    const double rn = r.norm();
    const Vector3d rh = r / rn;
    const double zeta = rh[2];
    const Vector3d dz(0.0, 0.0, 1.0);
    Matrix3d H = (15.0 * zeta * zeta - 3.0) * Matrix3d::Identity();
    H += 15.0 * rh * rh.transpose();
    H += -6.0 * dz * dz.transpose();
    H += 30.0 * zeta * (dz * rh.transpose() + rh * dz.transpose());
    H += -105.0 * zeta * zeta * rh * rh.transpose();
    return (C / std::pow(rn, 5)) * H;
}

Tensor3 j2_third(const Vector3d& r, double mu, double j2, double r_eq) {
    const double C = j2_C(mu, j2, r_eq);
    const double rn = r.norm();
    const Vector3d rh = r / rn;
    const double z = rh[2];                       // zeta = cos(theta)
    const double pref = C / std::pow(rn, 6);
    Tensor3 T;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) {
                const double diz = (i == 2), djz = (j == 2), dkz = (k == 2);
                const double S1 = (i == j) * rh[k] + (i == k) * rh[j] + (j == k) * rh[i];
                const double S2 = diz * djz * rh[k] + diz * dkz * rh[j] + djz * dkz * rh[i];
                const double S3 = diz * (j == k) + djz * (i == k) + dkz * (i == j);
                const double S4 = diz * rh[j] * rh[k] + djz * rh[i] * rh[k] + dkz * rh[i] * rh[j];
                T(i, j, k) = pref * (15.0 * (1.0 - 7.0 * z * z) * S1
                                     - 105.0 * (1.0 - 9.0 * z * z) * rh[i] * rh[j] * rh[k]
                                     + 30.0 * S2 + 30.0 * z * S3 - 210.0 * z * S4);
            }
    return T;
}

// ---- acceleration -------------------------------------------------------
static Vector3d kepler_accel(const Vector3d& u, double mu) {
    const double r = u.norm();
    return -mu * u / (r * r * r);
}

static Vector3d j2_accel(const Vector3d& r, double mu, double j2, double r_eq) {
    const double C = j2_C(mu, j2, r_eq);
    const double rn = r.norm();
    const double z = r[2];
    // a = -grad U_J2 ,  U_J2 = C(r^-3 - 3 z^2 r^-5)
    Vector3d g = -3.0 * C * r / std::pow(rn, 5) + 15.0 * C * z * z * r / std::pow(rn, 7);
    g[2] += -6.0 * C * z / std::pow(rn, 5);
    return -g;
}

Vector3d acceleration(const Vector3d& r, const PotentialModel& m) {
    Vector3d a = kepler_accel(r, m.central_gm);
    if (m.j2 != 0.0) a += j2_accel(r, m.central_gm, m.j2, m.r_eq);
    for (size_t l = 0; l < m.perturber_gm.size(); ++l)
        a += kepler_accel(r - m.perturber_pos[l], m.perturber_gm[l]);
    return a;
}

// ---- aggregates ---------------------------------------------------------
Matrix3d potential_hessian(const Vector3d& r, const PotentialModel& m) {
    Matrix3d H = kepler_hessian(r, m.central_gm);
    if (m.j2 != 0.0) H += j2_hessian(r, m.central_gm, m.j2, m.r_eq);
    for (size_t l = 0; l < m.perturber_gm.size(); ++l)
        H += kepler_hessian(r - m.perturber_pos[l], m.perturber_gm[l]);
    return H;
}

Tensor3 potential_third_derivative(const Vector3d& r, const PotentialModel& m) {
    Tensor3 T = kepler_third(r, m.central_gm);
    if (m.j2 != 0.0) T += j2_third(r, m.central_gm, m.j2, m.r_eq);
    for (size_t l = 0; l < m.perturber_gm.size(); ++l)
        T += kepler_third(r - m.perturber_pos[l], m.perturber_gm[l]);
    return T;
}

}  // namespace astdyn::propagation
