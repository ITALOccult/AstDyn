#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/time/TimeScale.hpp"
#include <iostream>

using namespace astdyn;

int main() {
    io::HorizonsClient horizons;
    auto mjd = time::calendar_to_mjd(2026, 3, 22, 12.58 / 24.0);
    auto t = time::EpochTDB::from_mjd(mjd);
    auto h = horizons.query_vectors("Vesta", t, "399"); 
    if (h) {
        Eigen::Vector3d rho = h->position.to_eigen_si();
        double r = rho.norm();
        double ra = std::atan2(rho.y(), rho.x()) * 180.0 / 3.14159265;
        double dec = std::asin(rho.z() / r) * 180.0 / 3.14159265;
        if (ra < 0) ra += 360.0;
        std::cout << "Vesta @ Earth Geometric RA: " << ra << " deg" << std::endl;
        std::cout << "Vesta @ Earth Geometric Dec: " << dec << " deg" << std::endl;
    }
    return 0;
}
