#include "astdyn/propagation/saba4_integrator.hpp"
#include "astdyn/core/Constants.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::constants;

int main() {
    std::cout << "--- SABA4 Reversibility Test (Two-Body) ---" << std::endl;

    // Apophis-like state (ICRF, AU-D)
    Eigen::VectorXd y0(6);
    y0 << -7.9651378683175167e-01, -5.1147558537755800e-01, -2.1018833794803629e-01,
           1.2113921469451630e-02, -1.1302603935421080e-02, -3.8947578481756361e-03;

    auto f = [](double t, const Eigen::VectorXd& y) {
        Eigen::VectorXd dy(6);
        double r3 = std::pow(y.head<3>().norm(), 3);
        dy.head<3>() = y.tail<3>();
        dy.tail<3>() = -(GMS / r3) * y.head<3>();
        return dy;
    };

    double h = 0.5; // step size in days
    SABA4Integrator saba(h, 1e-30, h, h);

    double E0 = 0.5 * y0.tail<3>().squaredNorm() - GMS / y0.head<3>().norm();
    
    double T_days = 50.0 * 365.25;
    std::cout << "Forward integration for 50 years (dt=" << h << ")..." << std::endl;
    Eigen::VectorXd yT = saba.integrate(f, y0, 0.0, T_days);

    double ET = 0.5 * yT.tail<3>().squaredNorm() - GMS / yT.head<3>().norm();
    double dE_fwd = std::abs((ET - E0) / E0);

    std::cout << "Reversing velocities..." << std::endl;
    Eigen::VectorXd y_rev = yT;
    y_rev.tail<3>() *= -1.0;

    std::cout << "Backward integration for 50 years..." << std::endl;
    Eigen::VectorXd y_back = saba.integrate(f, y_rev, 0.0, T_days);

    std::cout << "Re-reversing velocities..." << std::endl;
    y_back.tail<3>() *= -1.0;

    double E_back = 0.5 * y_back.tail<3>().squaredNorm() - GMS / y_back.head<3>().norm();
    double dE_total = std::abs((E_back - E0) / E0);

    double dr = (y_back.head<3>() - y0.head<3>()).norm();
    double eps_r = dr / y0.head<3>().norm();

    std::cout << std::scientific << std::setprecision(15);
    std::cout << "Initial r: " << y0.head<3>().transpose() << std::endl;
    std::cout << "Final r:   " << y_back.head<3>().transpose() << std::endl;
    std::cout << "Energy E0: " << E0 << std::endl;
    std::cout << "Rel Energy Error (fwd): " << dE_fwd << std::endl;
    std::cout << "Rel Energy Error (total): " << dE_total << std::endl;
    std::cout << "Error dr:  " << dr << " AU" << std::endl;
    std::cout << "Error eps_r: " << eps_r << std::endl;

    if (eps_r < 1e-11) {
        std::cout << "SUCCESS: SABA4 is reversible in Two-Body!" << std::endl;
    } else {
        std::cout << "FAILURE: SABA4 is NOT reversible (eps_r > 1e-11)" << std::endl;
    }

    return 0;
}
