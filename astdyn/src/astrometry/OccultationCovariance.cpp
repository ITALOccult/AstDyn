#include "astdyn/astrometry/OccultationCovariance.hpp"
#include <cmath>
#include <stdexcept>

namespace astdyn::astrometry {
namespace {

/// Unit line of sight and its length.
struct Direction { Vector3d rho_hat; double rho; };

Direction line_of_sight(const Vector3d& r, const Vector3d& r_obs) {
    const Vector3d d = r - r_obs;
    const double rho = d.norm();
    if (rho <= 0.0) throw std::runtime_error("OccultationCovariance: null line of sight");
    return { d / rho, rho };
}

} // namespace

TangentBasis tangent_basis(const Vector3d& s_hat, const Vector3d& pole) {
    const Vector3d s = s_hat.normalized();
    Vector3d e0 = pole.cross(s);
    if (e0.norm() < 1e-12) {
        throw std::runtime_error("OccultationCovariance: pole is parallel to the line of sight");
    }
    e0.normalize();
    const Vector3d e1 = s.cross(e0);
    return { e0, e1 };
}

ProjectionMaps projection_maps(const Vector3d& r, const TangentBasis& basis,
                               const Vector3d& r_obs) {
    const auto los = line_of_sight(r, r_obs);
    const Vector3d& u = los.rho_hat;
    const double rho = los.rho;

    ProjectionMaps pm;
    pm.A_h.setZero();
    pm.H_h[0].setZero();
    pm.H_h[1].setZero();

    // d2 rho_hat_i / dr_j dr_k, built once and contracted with each basis vector:
    //   (1/rho^2)(3 u_i u_j u_k - delta_ij u_k - delta_ik u_j - delta_jk u_i)
    for (int a = 0; a < 2; ++a) {
        const Vector3d& e = basis[a];

        // Eq. 4: (1/rho)(e - (e.u) u). Reduces to e/rho when e _|_ u, i.e. at the event.
        pm.A_h.block<1, 3>(a, 0) = ((e - (e.dot(u)) * u) / rho).transpose();

        // Eq. 7: contract e_i with d2 u_i / dr_j dr_k.
        const double eu = e.dot(u);
        Eigen::Matrix3d h;
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                h(j, k) = (3.0 * eu * u(j) * u(k)
                           - e(j) * u(k)
                           - e(k) * u(j)
                           - (j == k ? eu : 0.0)) / (rho * rho);
            }
        }
        pm.H_h[a].block<3, 3>(0, 0) = h;   // velocity blocks stay zero
    }
    return pm;
}

CompositeMap composite_map(const ProjectionMaps& pm, const Matrix6d& Phi,
                           const math::Tensor6& Psi) {
    CompositeMap cm;

    // Eq. 5
    cm.A = pm.A_h * Phi;

    // Eq. 13: B_ajk = A_h[a,p] Psi[p,j,k] + Phi^T H_h[a] Phi
    for (int a = 0; a < 2; ++a) {
        Matrix6d acc = Matrix6d::Zero();
        for (int p = 0; p < 6; ++p) {
            acc += pm.A_h(a, p) * Psi.slice(p);
        }
        acc += Phi.transpose() * pm.H_h[a] * Phi;
        cm.B[a] = acc;
    }
    return cm;
}

Moments moments(const CompositeMap& cm, const Matrix6d& C0) {
    Moments m;

    // Eq. 14: bias_a = 1/2 B_ajk C0_jk = 1/2 tr(B_a C0), since both are symmetric.
    for (int a = 0; a < 2; ++a) {
        m.bias(a) = 0.5 * (cm.B[a] * C0).trace();
    }

    // Eq. 15, linear part: A C0 A^T
    m.C_lin = cm.A * C0 * cm.A.transpose();

    // Eq. 15, second-order part: 1/2 B_ajk B_blm C0_jl C0_km = 1/2 tr(B_a C0 B_b C0)
    for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 2; ++b) {
            m.C_xi(a, b) = m.C_lin(a, b)
                         + 0.5 * (cm.B[a] * C0 * cm.B[b] * C0).trace();
        }
    }

    // Eq. 16: term_abc = A_aj A_bk B_clm C0_jl C0_km, then summed over the three
    // cyclic permutations. With v_a = C0 A_a^T this is v_a^T B_c v_b.
    std::array<Vector6d, 2> v;
    for (int a = 0; a < 2; ++a) v[a] = C0 * cm.A.row(a).transpose();

    std::array<std::array<Eigen::Vector2d, 2>, 2> term;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int c = 0; c < 2; ++c)
                term[a][b](c) = v[a].transpose() * cm.B[c] * v[b];

    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int c = 0; c < 2; ++c)
                m.kappa[a][b](c) = term[a][b](c) + term[b][c](a) + term[c][a](b);

    return m;
}

double nonlinearity_index(const CompositeMap& cm, const Matrix6d& C0,
                          const Eigen::Vector2d& n_hat) {
    // Project the maps onto n_hat.
    const Vector6d A_n = (n_hat.transpose() * cm.A).transpose();
    const Matrix6d B_n = n_hat(0) * cm.B[0] + n_hat(1) * cm.B[1];

    // RMS of the second-order term: sqrt(2 B_n:C0:B_n) = sqrt(2 tr(B_n C0 B_n C0))
    const double num = std::sqrt(2.0 * (B_n * C0 * B_n * C0).trace());
    // RMS of the first-order term: sqrt(A_n^T C0 A_n)
    const double den = std::sqrt(A_n.transpose() * C0 * A_n);
    if (den <= 0.0) throw std::runtime_error("OccultationCovariance: null first-order term");
    return 0.5 * num / den;
}

} // namespace astdyn::astrometry
