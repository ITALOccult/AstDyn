#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::astrometry;

int main() {
    io::HorizonsClient horizons;
    auto mjd = time::calendar_to_mjd(2026, 3, 22, 12.58265 / 24.0);
    auto t = time::EpochTDB::from_mjd(mjd);
    
    // 1. Vesta Apparent from Horizons (query observer)
    // RA/Dec Apparent usually includes light-time, aberration, but NOT refraction.
    // My HorizonsClient only has query_vectors. 
    // I will use vectors and manually apply light-time and aberration.
    
    // Star (from user)
    double star_ra = 22.59056205 * 15.0;
    double star_dec = -12.5368901;

    std::cout << "Target Star RA: " << star_ra << " deg, Dec: " << star_dec << " deg" << std::endl;

    auto h_geo = horizons.query_vectors("Vesta", t, "399"); // Geocentric
    if (h_geo) {
        Eigen::Vector3d rho = h_geo->position.to_eigen_si();
        double dist = rho.norm();
        double ra_geo = std::atan2(rho.y(), rho.x()) * 180.0 / 3.14159265;
        double dec_geo = std::asin(rho.z() / dist) * 180.0 / 3.14159265;
        if (ra_geo < 0) ra_geo += 360.0;
        
        std::cout << "Vesta Horizons Geometric RA: " << ra_geo << " deg, Dec: " << dec_geo << " deg" << std::endl;
        
        // Let's assume the star from the user is already Apparent.
        // Difference:
        double d_ra = (ra_geo - star_ra) * 3600.0 * std::cos(star_dec * 3.14159/180.0);
        double d_dec = (dec_geo - star_dec) * 3600.0;
        
        std::cout << "Geometric Separation: dRA*cosD=" << d_ra << "\" dDec=" << d_dec << "\"" << std::endl;
    }
    
    return 0;
}
