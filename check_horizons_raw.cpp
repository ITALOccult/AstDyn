#include "astdyn/io/HorizonsClient.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace std;

int main() {
    try {
        io::HorizonsClient horizons;
        time::EpochTDB t = time::EpochTDB::from_jd(2461171.1);
        
        cout << "JD: " << fixed << setprecision(6) << t.jd() << endl;
        
        auto elements = horizons.query_elements("234", t);
        if (elements) {
            cout << "Elements Semi-major axis: " << elements->a.to_au() << " AU" << endl;
            cout << "Elements eccentricity: " << elements->e << endl;
            cout << "Elements inclination: " << elements->i.to_deg() << " deg" << endl;
        } else {
            cout << "Failed to query elements." << endl;
        }

        // Check if there are other 234s
        // Usually, 234 is Barbara.
        
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }
    return 0;
}
