#include "astdyn_wrapper.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace ioccultcalc;

// Helper to convert RA (deg) to HMS
std::string degToHMS(double deg) {
    double h = deg / 15.0;
    int hh = static_cast<int>(h);
    double m = (h - hh) * 60.0;
    int mm = static_cast<int>(m);
    double s = (m - mm) * 60.0;
    
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << hh << " "
       << std::setw(2) << mm << " "
       << std::fixed << std::setprecision(2) << s;
    return ss.str();
}

// Helper to convert Dec (deg) to DMS
std::string degToDMS(double deg) {
    char sign = deg >= 0 ? '+' : '-';
    deg = std::abs(deg);
    int d = static_cast<int>(deg);
    double m = (deg - d) * 60.0;
    int mm = static_cast<int>(m);
    double s = (m - mm) * 60.0;
    
    std::stringstream ss;
    ss << sign << std::setfill('0') << std::setw(2) << d << " "
       << std::setw(2) << mm << " "
       << std::fixed << std::setprecision(1) << s;
    return ss.str();
}

// Helper to convert HMS string to degrees
double hmsToDeg(int h, int m, double s) {
    return (h + m/60.0 + s/3600.0) * 15.0;
}

// Helper to convert DMS string to degrees
double dmsToDeg(int d, int m, double s, bool negative) {
    double val = d + m/60.0 + s/3600.0;
    return negative ? -val : val;
}

struct TruthPoint {
    double mjd;
    std::string date_str;
    double ra_deg;
    double dec_deg;
};

int main() {
    std::cout << "================================================================================" << std::endl;
    std::cout << "  ASTEROID 249 (ILSE) HIGH-PRECISION VALIDATION (Jan 9-11 2026)" << std::endl;
    std::cout << "================================================================================" << std::endl;

    // 1. Initialize Wrapper
    AstDynWrapper wrapper;

    // 2. Set Osculating Elements (Fetched from JPL Horizons for Epoch 2461050.5 = MJD 61050.0? No wait, fetch said MJD 61050.5)
    // fetch_elements output: 2461051.0 = 2026-Jan-10 12:00 TDB -> MJD 61050.5 (if MJD=JD-2400000.5)
    double epoch_mjd = 61050.5;
    
    // Elements from fetch_elements.py output:
    // A = 3.558083341572964E+08 km
    double au_km = 149597870.7;
    double a_au = 3.558083341572964e8 / au_km; 
    
    // EC= 2.173906694003822E-01
    double e = 0.2173906694003822;
    
    // IN= 9.621856552411442E+00
    double i_deg = 9.621856552411442;
    
    // OM= 3.346468923829066E+02
    double Omega_deg = 334.6468923829066;
    
    // W = 4.237096059969402E+01
    double omega_deg = 42.37096059969402;
    
    // MA= 7.165173885188509E+01
    double M_deg = 71.65173885188509;

    // Convert angles to radians
    double deg2rad = M_PI / 180.0;
    
    // Initialize wrapper with precise Keplerian elements
    // Note: setKeplerianElements expects standard Keplerian (Ecliptic J2000 usually)
    // JPL Horizons elements for 249 are Ecliptic J2000.
    wrapper.setKeplerianElements(
        a_au, e, i_deg * deg2rad,
        Omega_deg * deg2rad, omega_deg * deg2rad, M_deg * deg2rad,
        epoch_mjd, "249 Ilse",
        astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC
    );

    std::cout << "Initialized Elements at MJD " << epoch_mjd << " (2026-Jan-10 12:00 TDB)" << std::endl;
    std::cout << "  a: " << a_au << " AU" << std::endl;
    std::cout << "  e: " << e << std::endl;
    std::cout << "  i: " << i_deg << " deg" << std::endl;

    // 3. Truth Data (Jan 9, 10, 11)
    std::vector<TruthPoint> truth_points;
    
    // Jan 09 00:00 (MJD 61049.0): 08 08 43.06 +31 18 55.5
    truth_points.push_back({61049.0, "2026-Jan-09", hmsToDeg(8, 8, 43.06), dmsToDeg(31, 18, 55.5, false)});
    
    // Jan 10 00:00 (MJD 61050.0): 08 07 25.53 +31 19 21.6
    truth_points.push_back({61050.0, "2026-Jan-10", hmsToDeg(8, 7, 25.53), dmsToDeg(31, 19, 21.6, false)});
    
    // Jan 11 00:00 (MJD 61051.0): 08 06 07.45 +31 19 38.2
    truth_points.push_back({61051.0, "2026-Jan-11", hmsToDeg(8, 6, 7.45), dmsToDeg(31, 19, 38.2, false)});

    // 4. Verification Loop
    std::cout << "\nValidation Results:" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << " Date        | Ref MJD | RA (Calc)     | RA (Truth)    | dRA('')| Dec(Calc)     | Dec(Truth)    | dDec('')" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;

    double max_err_ra = 0.0;
    double max_err_dec = 0.0;

    for (const auto& tp : truth_points) {
        // Calculate observation
        auto obs = wrapper.calculateObservation(tp.mjd);
        
        double d_ra = (obs.ra_deg - tp.ra_deg) * 3600.0;
        double d_dec = (obs.dec_deg - tp.dec_deg) * 3600.0;
        
        if (std::abs(d_ra) > max_err_ra) max_err_ra = std::abs(d_ra);
        if (std::abs(d_dec) > max_err_dec) max_err_dec = std::abs(d_dec);

        std::cout << tp.date_str << " | " 
                  << std::fixed << std::setprecision(1) << tp.mjd << " | "
                  << degToHMS(obs.ra_deg) << " | " << degToHMS(tp.ra_deg) << " | " 
                  << std::setw(6) << std::setprecision(3) << d_ra << " | "
                  << degToDMS(obs.dec_deg) << " | " << degToDMS(tp.dec_deg) << " | " 
                  << std::setw(6) << std::setprecision(3) << d_dec << std::endl;
    }
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    
    std::cout << "\nMAX ERROR -> RA: " << max_err_ra << " arcsec, Dec: " << max_err_dec << " arcsec" << std::endl;
    
    if (max_err_ra < 1.0 && max_err_dec < 1.0) {
        std::cout << "\n[SUCCESS] Precision < 1 arcsec verified!" << std::endl;
        return 0;
    } else {
        std::cout << "\n[WARNING] Precision check failed threshold (< 1 arcsec)." << std::endl;
        return 1;
    }
}
