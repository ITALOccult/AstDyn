#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <iostream>
#include <iomanip>
using namespace astdyn;
int main(){
    // Greenwich in ITRF: (R, 0, 0)
    Eigen::Vector3d g_itrf(6378137.0, 0, 0);
    for (double mjd : {51544.5, 51544.5+0.25, 51544.5+0.5}) {
        auto t = time::EpochUTC::from_mjd(mjd);
        double era = coordinates::ReferenceFrame::computeERA(t, 0.0);
        auto R = coordinates::ReferenceFrame::itrf_to_j2000_simple(t);
        Eigen::Vector3d g = R * g_itrf;
        std::cout << std::fixed << std::setprecision(4)
                  << "MJD " << mjd << "  ERA=" << era*180/M_PI << " deg\n"
                  << "   calcolato  y = " << g.y()/1000 << " km\n"
                  << "   atteso     y = " << 6378.137*std::sin(era) << " km\n";
    }
}
