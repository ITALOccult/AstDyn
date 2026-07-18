/**
 * @file test_force_field_gradient.cpp
 * @brief Valida ForceField::acceleration_gradient / _second_gradient contro
 *        le differenze finite di total_acceleration (giudice indipendente
 *        dalle formule analitiche di derivata).
 */
#include <gtest/gtest.h>
#include "astdyn/propagation/ForceField.hpp"
#include "astdyn/propagation/PropagatorSettings.hpp"

using namespace astdyn;
using namespace astdyn::propagation;

// Settings kepler+J2 solare, N-body spento -> nessuna effemeride necessaria.
static PropagatorSettings kepler_j2_settings() {
    PropagatorSettings s;
    s.include_planets   = false;
    s.include_asteroids = false;
    s.include_relativity= false;
    s.include_earth_j2  = false;
    s.include_sun_j2    = true;
    s.baricentric_integration = false;   // eliocentrico
    return s;
}

// Gradiente numerico: dA_i/dr_j via differenze centrate su total_acceleration.
static Eigen::Matrix3d fd_gradient(const ForceField& ff, time::EpochTDB t,
                                   const Eigen::Vector3d& r, double h) {
    Eigen::Matrix3d G;
    const Eigen::Vector3d v = Eigen::Vector3d::Zero();
    for (int j = 0; j < 3; ++j) {
        Eigen::Vector3d rp = r, rm = r; rp[j] += h; rm[j] -= h;
        Eigen::Vector3d d = (ff.total_acceleration(t, rp, v)
                           - ff.total_acceleration(t, rm, v)) / (2*h);
        for (int i = 0; i < 3; ++i) G(i, j) = d[i];
    }
    return G;
}

TEST(ForceFieldGradient, AccelGradientMatchesFiniteDiff) {
    ForceField ff(kepler_j2_settings());
    time::EpochTDB t = time::EpochTDB::from_mjd(60000.0);
    const double h = 1e-7;
    const Eigen::Vector3d pts[] = {
        {1.6, 0.3, 0.2}, {2.1, -1.0, 0.7}, {0.9, 0.9, -0.5}};
    for (const auto& r : pts) {
        Eigen::Matrix3d Gana = ff.acceleration_gradient(t, r, Eigen::Vector3d::Zero());
        Eigen::Matrix3d Gnum = fd_gradient(ff, t, r, h);
        double rel = (Gana - Gnum).norm() / Gnum.norm();
        EXPECT_LT(rel, 1e-5) << "punto (" << r.transpose() << ")  rel=" << rel;
    }
}

TEST(ForceFieldGradient, SecondGradientMatchesFiniteDiffOfGradient) {
    ForceField ff(kepler_j2_settings());
    time::EpochTDB t = time::EpochTDB::from_mjd(60000.0);
    const double h = 1e-6;
    const Eigen::Vector3d r(1.6, 0.3, 0.2);
    // T_ijk = dG_ij/dr_k  (derivata numerica del gradiente analitico)
    math::Tensor3 Tana = ff.acceleration_second_gradient(t, r, Eigen::Vector3d::Zero());
    double num = 0, den = 0;
    for (int k = 0; k < 3; ++k) {
        Eigen::Vector3d rp = r, rm = r; rp[k] += h; rm[k] -= h;
        Eigen::Matrix3d dG =
            (ff.acceleration_gradient(t, rp, Eigen::Vector3d::Zero())
           - ff.acceleration_gradient(t, rm, Eigen::Vector3d::Zero())) / (2*h);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double d = Tana(i,j,k) - dG(i,j);
                num += d*d; den += dG(i,j)*dG(i,j);
            }
    }
    EXPECT_LT(std::sqrt(num/den), 1e-4);
}
