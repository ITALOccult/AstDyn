#include "astdyn_wrapper.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace ioccultcalc;

// Dati di verità per il test (copiati da test_34713_precision.cpp, Case A result)
// MJD 61050.0
// RA:  85.96338 deg (calcolato da Case A)
// Dec: 26.51629 deg (calcolato da Case A)
// Vogliamo verificare che usando setMeanEquinoctialElements otteniamo lo stesso risultato.

int main() {
    std::cout << "================================================================================" << std::endl;
    std::cout << "  ALLNUM.CAT DB SIMULATION TEST (Mean Elements -> High Precision)" << std::endl;
    std::cout << "================================================================================" << std::endl;

    AstDynWrapper wrapper;

    // --- SIMULATED DB RECORD (Asteroid 34713) ---
    // Estratto manualmente da astdyn/data/34713.eq1 (Mean Equinoctial Elements)
    // EPOCH_MJD = 58016.5
    // a = 2.377236879199628
    // h = 0.1601006093284051
    // k = 0.1477430150931123
    // p = -0.01633534579973804
    // q = -0.08253579899182363
    // lambda = 24.32289433439971 (Mean Longitude in degrees? No, radians usually in code, check file format)
    // .eq1 FILE: 
    // ML=  24.32289433439971 (Mean Longitude is usually in degrees in eq1 file, need to interpret correctly)
    
    // Let's re-read how parse_eq1 works or check the values used.
    // parse_eq1 uses OrbFitEQ1Parser.
    // EQUINOCTIAL ELEMENTS in .eq1 are:
    // a (AU)
    // h = e * sin(pomega)
    // k = e * cos(pomega)
    // p = tan(i/2) * sin(Omega)
    // q = tan(i/2) * cos(Omega)
    // lambda = M + pomega (Mean Longitude)
    
    // Values from 34713.eq1:
    // A=  2.377236879199628
    // H=  0.1601006093284051
    // K=  0.1477430150931123
    // P= -0.01633534579973804
    // Q= -0.08253579899182363
    // ML= 24.32289433439971  (Degrees? Parser usually handles rad/deg conversion. Let's assume input to wrapper is RAD for angles)
    
    // Values from 34713.eq1:
    // EQU   2.8041167798722424E+00  -0.018233692083002   0.172196556501974   -0.065799001000219    0.023094080542415  63.7347677138408
    // MJD     61000.000000000 TDT
    
    double epoch_mjd = 61000.0;
    double a = 2.8041167798722424;
    double h = -0.018233692083002;
    double k = 0.172196556501974;
    double p = -0.065799001000219;
    double q = 0.023094080542415;
    double ml_deg = 63.7347677138408; // Degrees in .eq1
    double ml_rad = ml_deg * M_PI / 180.0;

    std::cout << "Setting Mean Equinoctial Elements (Simulating DB fetch)..." << std::endl;
    // Input to setMeanEquinoctialElements expects Radians for lambda? 
    // Let's check setEquinoctialElements: expects lambda [rad].
    // So we pass ml_rad.
    
    wrapper.setMeanEquinoctialElements(a, h, k, p, q, ml_rad, epoch_mjd, "34713 (DB)");
    
    // Verify initialization
    if (!wrapper.isInitialized()) {
        std::cerr << "ERROR: Wrapper failed to initialize!" << std::endl;
        return 1;
    }
    
    // --- VERIFICATION (Against JPL Truth 61050.0) ---
    // From test_34713_precision.cpp Case A:
    // RA:  85.9633827299 (Diff: 0.627 arcsec)
    // Dec: 26.5162919830 (Diff: 0.151 arcsec)
    
    double target_mjd = 61050.0;
    std::cout << "Propagating to MJD " << target_mjd << "..." << std::endl;
    
    auto result = wrapper.calculateObservation(target_mjd);
    
    std::cout << "Result:" << std::endl;
    std::cout << "  RA:  " << std::fixed << std::setprecision(8) << result.ra_deg << " deg" << std::endl;
    std::cout << "  Dec: " << std::fixed << std::setprecision(8) << result.dec_deg << " deg" << std::endl;
    
    // Validate against Expected (Case A result)
    double expected_ra = 85.9633827; 
    double expected_dec = 26.5162919;
    
    double diff_ra = std::abs(result.ra_deg - expected_ra) * 3600.0;
    double diff_dec = std::abs(result.dec_deg - expected_dec) * 3600.0;
    
    std::cout << "Difference from Baseline (File Loading):" << std::endl;
    std::cout << "  dRA:  " << diff_ra << " arcsec" << std::endl;
    std::cout << "  dDec: " << diff_dec << " arcsec" << std::endl;
    
    if (diff_ra < 0.001 && diff_dec < 0.001) {
        std::cout << "[SUCCESS] DB Simulation matches File Loading perfectly." << std::endl;
        return 0;
    } else {
        std::cout << "[FAILURE] DB Simulation result deviates by > 0.001 arcsec." << std::endl;
        return 1;
    }
}
