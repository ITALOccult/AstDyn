#include "astdyn/AstDyn.hpp"
#include <curl/curl.h>
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace std;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.ephemeris_type = EphemerisType::DE441;
        
        double jd_utc = 2461164.5;
        auto t_utc = time::EpochTDB::from_jd(jd_utc); 
        
        io::HorizonsClient horizons;
        
        // 1. Query HELIOCENTRIC ELEMENTS from JPL
        auto kep_opt = horizons.query_elements("234", t_utc);
        if (!kep_opt) return 1;
        auto elements = kep_opt.value();
        
        cout << fixed << setprecision(10);
        cout << "Elements parsed: A=" << elements.a.to_au() << " E=" << elements.e << " I=" << elements.i.to_deg() << endl;
        cout << "OM=" << elements.node.to_deg() << " W=" << elements.omega.to_deg() << " MA=" << elements.M.to_deg() << endl;

        // 2. Convert to Cartesian
        auto cart_helio_ecl = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(elements);
        auto p_local = cart_helio_ecl.position.to_eigen_si();
        cout << "Local Cartesian (Helio Ecl, m): " << p_local.transpose() << endl;
        
        // 3. Query JPL VECTORS (Heliocentric Ecliptic)
        string url = "https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND='234'&OBJ_DATA='NO'&MAKE_EPHEM='YES'&EPHEM_TYPE='VECTORS'&CENTER='@10'&START_TIME='JD2461164.5'&STOP_TIME='JD2461164.6'&STEP_SIZE='1d'&REF_PLANE='ECLIPTIC'&REF_SYSTEM='J2000'&CSV_FORMAT='NO'&VEC_TABLE='2'";
        
        CURL* curl = curl_easy_init();
        string buffer;
        if (curl) {
            curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
            auto write_cb = [](void* c, size_t s, size_t n, void* u) {
                ((string*)u)->append((char*)c, s * n); return s * n;
            };
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, +write_cb);
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, &buffer);
            curl_easy_perform(curl);
            curl_easy_cleanup(curl);
        }
        
        cout << "\nJPL API Response for Heliocentric Ecliptic:\n";
        size_t start = buffer.find("$$SOE");
        if (start != string::npos) {
            size_t end = buffer.find("$$EOE");
            if (end != string::npos) cout << buffer.substr(start, end - start + 5) << endl;
            else cout << buffer.substr(start) << endl;
        }

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
