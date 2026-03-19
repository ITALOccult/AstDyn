#include "astdyn/AstDyn.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace astdyn;
using namespace std;

int main() {
    try {
        io::HorizonsClient horizons;
        cout << fixed << setprecision(6);
        cout << "UTC_Time, JD, RA_deg, Dec_deg, Dist_AU" << endl;
        
        for (int d = 1; d <= 31; ++d) {
            double jd = 2461161.5 + d; // May 1 is JD 2461161.5
            time::EpochTDB t = time::EpochTDB::from_jd(jd);
            
            // Query Astrometric RA/Dec from Geocenter (500)
            auto obs = horizons.query_observation("234", t, "500");
            if (obs) {
                cout << "2026-05-" << (d < 10 ? "0" : "") << d << " 00:00:00, "
                     << jd << ", " << obs->ra.to_deg() << ", " << obs->dec.to_deg() << ", " << obs->distance.to_au() << endl;
            }
        }
    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
