#include <astdyn/AstDynEngine.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/time/TimeScale.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

using namespace astdyn;

/**
 * @brief Multi-body Dynamics for Haumea System
 * States: [H_pos, H_vel, Hi_pos, Hi_vel, Na_pos, Na_vel] (18 elements)
 * Origin: Heliocentric (Sun at 0,0,0)
 */
struct HaumeaSystemDynamics {
    double gm_sun = constants::GM_SUN;
    double gm_h = 267.4;    // Haumea (km^3/s^2)
    double gm_hi = 1.2;     // Hi'iaka
    double gm_na = 0.12;    // Namaka

    Eigen::VectorXd operator()(double t, const Eigen::VectorXd& y) const {
        Eigen::VectorXd dy(18);
        
        Eigen::Vector3d r_h  = y.segment<3>(0);
        Eigen::Vector3d v_h  = y.segment<3>(3);
        Eigen::Vector3d r_hi = y.segment<3>(6);
        Eigen::Vector3d v_hi = y.segment<3>(9);
        Eigen::Vector3d r_na = y.segment<3>(12);
        Eigen::Vector3d v_na = y.segment<3>(15);

        // Accelerations from Sun
        Eigen::Vector3d a_h_sun  = -gm_sun * r_h  / std::pow(r_h.norm(), 3);
        Eigen::Vector3d a_hi_sun = -gm_sun * r_hi / std::pow(r_hi.norm(), 3);
        Eigen::Vector3d a_na_sun = -gm_sun * r_na / std::pow(r_na.norm(), 3);

        // Mutual attractions
        auto acc_mutual = [&](const Eigen::Vector3d& r1, const Eigen::Vector3d& r2, double gm2) {
            Eigen::Vector3d dist = r1 - r2;
            return -gm2 * dist / std::pow(dist.norm(), 3);
        };

        Eigen::Vector3d a_h_mutual  = acc_mutual(r_h, r_hi, gm_hi) + acc_mutual(r_h, r_na, gm_na);
        Eigen::Vector3d a_hi_mutual = acc_mutual(r_hi, r_h, gm_h)  + acc_mutual(r_hi, r_na, gm_na);
        Eigen::Vector3d a_na_mutual = acc_mutual(r_na, r_h, gm_h)  + acc_mutual(r_na, r_hi, gm_hi);

        dy.segment<3>(0) = v_h;
        dy.segment<3>(3) = a_h_sun + a_h_mutual;
        dy.segment<3>(6) = v_hi;
        dy.segment<3>(9) = a_hi_sun + a_hi_mutual;
        dy.segment<3>(12) = v_na;
        dy.segment<3>(15) = a_na_sun + a_na_mutual;

        return dy;
    }
};

int main() {
    std::cout << "=== Haumea Multi-Body System Evolution ===\n";

    // Initial States (JD 2461164.5, km and km/s, Ecliptic J2000)
    Eigen::VectorXd y0(18);
    
    // Haumea Bary (Heliocentric)
    y0(0) = -5.513013154703953e+09; y0(1) = -3.563873466876123e+09; y0(2) = +3.520538626605666e+09;
    y0(3) = +2.411700979939033e+00; y0(4) = -3.024644217319804e+00; y0(5) = -2.449900086844126e-01;
    
    // Hi'iaka (Barycentric -> Heliocentric)
    Eigen::Vector3d hi_rel(+3.230691247696847e+04, +2.899365818987711e+04, +1.578694964659445e+04);
    Eigen::Vector3d hi_v_rel(+4.835737731380611e-02, -2.198936577535662e-02, -5.659653697591440e-02);
    y0.segment<3>(6) = y0.segment<3>(0) + hi_rel;
    y0.segment<3>(9) = y0.segment<3>(3) + hi_v_rel;

    // Namaka (Barycentric -> Heliocentric)
    Eigen::Vector3d na_rel(+2.117923575673045e+04, +5.998545766923913e+03, -1.132606650049148e+04);
    Eigen::Vector3d na_v_rel(-2.695122300908965e-02, -8.481368036922111e-02, -5.418276147615637e-02);
    y0.segment<3>(12) = y0.segment<3>(0) + na_rel;
    y0.segment<3>(15) = y0.segment<3>(3) + na_v_rel;

    // Integration setup
    double t0 = 0;
    double tf = 50.0 * 86400.0; // 50 days
    double dt_out = 0.5 * 86400.0; // Sample every 12 hours
    
    HaumeaSystemDynamics dynamics;
    propagation::RKF78Integrator integrator(600.0, 1e-13);
    
    std::vector<double> t_steps;
    std::vector<Eigen::VectorXd> y_steps;
    
    std::cout << "Propagating system for 50 days (RKF78)...\n";
    for (double t = t0; t <= tf; t += dt_out) {
        t_steps.push_back(t);
        if (t == t0) {
            y_steps.push_back(y0);
        } else {
            y_steps.push_back(integrator.integrate(dynamics, y_steps.back(), t - dt_out, t));
        }
    }

    // Export relative positions
    std::ofstream out("haumea_system_evolution.csv");
    out << "time_days,hi_x_rel,hi_y_rel,hi_z_rel,hi_dist,na_x_rel,na_y_rel,na_z_rel,na_dist\n";
    for (size_t i = 0; i < t_steps.size(); ++i) {
        double days = t_steps[i] / 86400.0;
        Eigen::Vector3d r_h = y_steps[i].segment<3>(0);
        Eigen::Vector3d r_hi = y_steps[i].segment<3>(6);
        Eigen::Vector3d r_na = y_steps[i].segment<3>(12);
        
        Eigen::Vector3d d_hi = r_hi - r_h;
        Eigen::Vector3d d_na = r_na - r_h;
        
        out << days << "," 
            << d_hi.x() << "," << d_hi.y() << "," << d_hi.z() << "," << d_hi.norm() << ","
            << d_na.x() << "," << d_na.y() << "," << d_na.z() << "," << d_na.norm() << "\n";
    }
    out.close();

    std::cout << "Results saved to haumea_system_evolution.csv\n";
    
    // Summary
    std::cout << "\nFinal state (Day 50):\n";
    Eigen::Vector3d d_hi_f = y_steps.back().segment<3>(6) - y_steps.back().segment<3>(0);
    Eigen::Vector3d d_na_f = y_steps.back().segment<3>(12) - y_steps.back().segment<3>(0);
    std::cout << "Hi'iaka Distance: " << d_hi_f.norm() << " km\n";
    std::cout << "Namaka Distance: " << d_na_f.norm() << " km\n";

    return 0;
}
