#include "astdyn/core/physics_types.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "astdyn/astrometry/sky_types.hpp"
#include "src/core/frame_tags.hpp"
#include <iostream>

using namespace astdyn;
using namespace astdyn::physics;
using namespace astdyn::time;
using namespace astdyn::math;
using namespace astdyn::astrometry;

// Forward declare frames if needed 
namespace astdyn::core { struct GCRF; struct ECLIPJ2000; }

int main() {
    std::cout << "--- AstDyn 3.0 Type Safety Test ---\n";

    // 1. Physics Units
    Distance d1 = Distance::from_au(1.0);
    Velocity v1 = Velocity::from_km_s(30.0);
    GravitationalParameter gm = GravitationalParameter::sun();
    
    std::cout << "1 AU in meters: " << d1.to_m() << "\n";
    std::cout << "GM_SUN in m3/s2: " << gm.to_m3_s2() << "\n";
    // d1 + v1;  // <--- WOULD NOT COMPILE!

    // 2. Epochs and TimeDuration
    EpochTDB t1 = EpochTDB::from_jd(2451545.0); // J2000
    EpochTDB t2 = t1 + TimeDuration::from_days(10.0);
    TimeDuration delta = t2 - t1;
    
    std::cout << "Delta time in seconds: " << delta.to_seconds() << "\n";
    
    // EpochUTC t_utc = EpochUTC::from_jd(2451545.0);
    // t_utc - t1; // <--- WOULD NOT COMPILE! (different scales)

    // 3. Frame-bound Vectors and Rotation
    Vector3<core::ECLIPJ2000, Distance> r_ecl = Vector3<core::ECLIPJ2000, Distance>::from_si(1.5e11, 0, 0); // 1 AU x-axis
    
    // Fake rotation matrix (identity for test)
    Eigen::Matrix3d mat = Eigen::Matrix3d::Identity();
    RotationMatrix<core::ECLIPJ2000, core::GCRF> R_ecl_to_eq = RotationMatrix<core::ECLIPJ2000, core::GCRF>::from_eigen(mat);
    
    // Valid rotation
    Vector3<core::GCRF, Distance> r_eq = R_ecl_to_eq * r_ecl;
    
    // R_ecl_to_eq * r_eq; // <--- WOULD NOT COMPILE! (matrix expects ECLIPJ2000 input)

    std::cout << "Vector magnitude (m): " << r_eq.norm().to_m() << "\n";

    // 4. Astrometry
    // SkyCoord<core::GCRF> obs = SkyCoord<core::GCRF>::from_vector(r_ecl); // <--- WOULD NOT COMPILE! (needs GCRF)
    
    auto obs = SkyCoord<core::GCRF>::from_vector(r_eq);
    std::cout << "Right Ascension: " << obs.ra().to_deg() << " deg\n";
    
    std::cout << "--- All tests passed! ---\n";
    return 0;
}
