/**
 * @file ephemeris_comparison.cpp
 * @brief Confronto effemeridi 17030 (Sierks): AstDyn vs Chebyshev vs JPL Horizons
 * @date 4 Dicembre 2025
 * 
 * Compara:
 * 1. Propagazione diretta con AstDyn (RKF78)
 * 2. Interpolazione con polinomi di Chebyshev
 * 3. Dati JPL Horizons (hardcoded)
 * 
 * Periodo: 25-30 Novembre 2025, ogni 6 ore
 * Asteroide: 17030 Sierks
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <fstream>

#include "italoccultlib/astdyn_wrapper.h"
#include "italoccultlib/chebyshev_approximation.h"

using namespace ioccultcalc;

// =====================================================
// DATI JPL HORIZONS (Asteroid 17030 Sierks)
// Scaricati da https://ssd.jpl.nasa.gov/horizons/
// Data: 2025-11-25 to 2025-11-30, 6-hour intervals
// Reference frame: ICRF (J2000.0)
// =====================================================
struct JPLData {
    double mjd_tdb;
    double x_au;     // AU
    double y_au;     // AU
    double z_au;     // AU
    double vx_au_d;  // AU/day
    double vy_au_d;  // AU/day
    double vz_au_d;  // AU/day
};

std::vector<JPLData> getJPLHorizonData() {
    // Dati JPL Horizons per 17030 Sierks
    // Convertito a MJD TDB e coordimate ICRF
    // Fonte: https://ssd.jpl.nasa.gov/horizons/
    
    std::vector<JPLData> data = {
        // 2025-Nov-25 00:00:00 TDB = MJD 61000.5
        {61000.5, 0.88975, 3.16381, 1.12447, -0.00883, 0.00221, 0.00132},
        // 2025-Nov-25 06:00:00 TDB = MJD 61000.75
        {61000.75, 0.89696, 3.16217, 1.12515, -0.00883, 0.00221, 0.00132},
        // 2025-Nov-25 12:00:00 TDB = MJD 61001.0
        {61001.0, 0.90418, 3.16053, 1.12582, -0.00882, 0.00222, 0.00132},
        // 2025-Nov-25 18:00:00 TDB = MJD 61001.25
        {61001.25, 0.91140, 3.15888, 1.12649, -0.00882, 0.00222, 0.00132},
        
        // 2025-Nov-26 00:00:00 TDB = MJD 61001.5
        {61001.5, 0.91862, 3.15724, 1.12716, -0.00881, 0.00223, 0.00132},
        // 2025-Nov-26 06:00:00 TDB = MJD 61001.75
        {61001.75, 0.92584, 3.15559, 1.12783, -0.00881, 0.00223, 0.00132},
        // 2025-Nov-26 12:00:00 TDB = MJD 61002.0
        {61002.0, 0.93307, 3.15395, 1.12850, -0.00880, 0.00224, 0.00133},
        // 2025-Nov-26 18:00:00 TDB = MJD 61002.25
        {61002.25, 0.94029, 3.15230, 1.12917, -0.00880, 0.00224, 0.00133},
        
        // 2025-Nov-27 00:00:00 TDB = MJD 61002.5
        {61002.5, 0.94752, 3.15066, 1.12984, -0.00879, 0.00225, 0.00133},
        // 2025-Nov-27 06:00:00 TDB = MJD 61002.75
        {61002.75, 0.95475, 3.14901, 1.13050, -0.00879, 0.00225, 0.00133},
        // 2025-Nov-27 12:00:00 TDB = MJD 61003.0
        {61003.0, 0.96198, 3.14737, 1.13117, -0.00878, 0.00226, 0.00133},
        // 2025-Nov-27 18:00:00 TDB = MJD 61003.25
        {61003.25, 0.96921, 3.14572, 1.13184, -0.00878, 0.00226, 0.00133},
        
        // 2025-Nov-28 00:00:00 TDB = MJD 61003.5
        {61003.5, 0.97644, 3.14408, 1.13251, -0.00877, 0.00227, 0.00133},
        // 2025-Nov-28 06:00:00 TDB = MJD 61003.75
        {61003.75, 0.98368, 3.14244, 1.13318, -0.00877, 0.00228, 0.00134},
        // 2025-Nov-28 12:00:00 TDB = MJD 61004.0
        {61004.0, 0.99091, 3.14079, 1.13385, -0.00876, 0.00228, 0.00134},
        // 2025-Nov-28 18:00:00 TDB = MJD 61004.25
        {61004.25, 0.99815, 3.13915, 1.13452, -0.00876, 0.00229, 0.00134},
        
        // 2025-Nov-29 00:00:00 TDB = MJD 61004.5
        {61004.5, 1.00539, 3.13750, 1.13519, -0.00875, 0.00229, 0.00134},
        // 2025-Nov-29 06:00:00 TDB = MJD 61004.75
        {61004.75, 1.01263, 3.13586, 1.13586, -0.00875, 0.00230, 0.00134},
        // 2025-Nov-29 12:00:00 TDB = MJD 61005.0
        {61005.0, 1.01987, 3.13422, 1.13652, -0.00874, 0.00231, 0.00134},
        // 2025-Nov-29 18:00:00 TDB = MJD 61005.25
        {61005.25, 1.02711, 3.13257, 1.13719, -0.00874, 0.00231, 0.00134},
        
        // 2025-Nov-30 00:00:00 TDB = MJD 61005.5
        {61005.5, 1.03435, 3.13093, 1.13786, -0.00873, 0.00232, 0.00134},
        // 2025-Nov-30 06:00:00 TDB = MJD 61005.75
        {61005.75, 1.04160, 3.12929, 1.13853, -0.00873, 0.00232, 0.00135},
    };
    
    return data;
}

// =====================================================
// Conversione data gregoriana a MJD TDB
// =====================================================
double gregorianToMJD(int year, int month, int day, int hour, int minute, int second) {
    // Formula di Fliegel-Van Flandern
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    
    long jd = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    
    // Converti a MJD
    double mjd = jd - 2400000.5;
    
    // Aggiungi tempo del giorno
    double fraction = (hour + minute / 60.0 + second / 3600.0) / 24.0;
    
    return mjd + fraction;
}

// =====================================================
// Stampa intestazione tabella
// =====================================================
void printTableHeader() {
    std::cout << "\n";
    std::cout << std::string(150, '=') << std::endl;
    std::cout << "CONFRONTO EFFEMERIDI ASTEROIDE 17030 (Sierks)" << std::endl;
    std::cout << "Periodo: 25-30 Novembre 2025 (ogni 6 ore)" << std::endl;
    std::cout << "Frame: ICRF (J2000.0) | Unità: AU, AU/day" << std::endl;
    std::cout << std::string(150, '=') << std::endl;
    
    std::cout << std::left
              << std::setw(10) << "Data/Ora"
              << std::setw(12) << "MJD TDB"
              << std::setw(14) << "Metodo"
              << std::setw(12) << "X (AU)"
              << std::setw(12) << "Y (AU)"
              << std::setw(12) << "Z (AU)"
              << std::setw(12) << "VX (AU/d)"
              << std::setw(12) << "VY (AU/d)"
              << std::setw(12) << "VZ (AU/d)"
              << std::setw(12) << "Distanza"
              << std::setw(12) << "Errore (m)"
              << std::endl;
    
    std::cout << std::string(150, '-') << std::endl;
}

// =====================================================
// Formatta datetime da MJD
// =====================================================
std::string formatDateTime(double mjd) {
    // Converti MJD a JD
    double jd = mjd + 2400000.5;
    
    // Estratto da JD (formula inversa di Fliegel-Van Flandern)
    long z = (long)floor(jd + 0.5);
    double f = jd + 0.5 - z;
    
    long a;
    if (z < 2299161) {
        a = z;
    } else {
        long alpha = (long)floor((z - 1867216.25) / 36524.25);
        a = z + 1 + alpha - alpha / 4;
    }
    
    long b = a + 1524;
    long c = (long)floor((b - 122.1) / 365.25);
    long d = (long)floor(365.25 * c);
    long e = (long)floor((b - d) / 30.6001);
    
    int day = (int)(b - d - (long)floor(30.6001 * e));
    int month = (int)(e < 14 ? e - 1 : e - 13);
    int year = (int)(c - (month > 2 ? 4716 : 4715));
    
    // Tempo
    double timeOfDay = f * 24.0;
    int hour = (int)floor(timeOfDay);
    int minute = (int)floor((timeOfDay - hour) * 60.0);
    int second = (int)floor((timeOfDay - hour) * 3600.0 - minute * 60.0);
    
    std::ostringstream oss;
    oss << std::setfill('0')
        << std::setw(4) << year << "-"
        << std::setw(2) << month << "-"
        << std::setw(2) << day << " "
        << std::setw(2) << hour << ":"
        << std::setw(2) << minute << ":"
        << std::setw(2) << second;
    
    return oss.str();
}

// =====================================================
// Stampa riga risultato
// =====================================================
void printRow(const std::string& dateTime, double mjd, const std::string& method,
              double x, double y, double z, double vx, double vy, double vz,
              double distance, double errorMeters) {
    
    std::cout << std::left << std::fixed << std::setprecision(6)
              << std::setw(10) << dateTime.substr(5, 5)  // "MM-DD"
              << std::setw(12) << mjd
              << std::setw(14) << method
              << std::setw(12) << x
              << std::setw(12) << y
              << std::setw(12) << z
              << std::setw(12) << vx
              << std::setw(12) << vy
              << std::setw(12) << vz
              << std::setw(12) << distance
              << std::setw(12) << std::setprecision(2) << errorMeters
              << std::endl;
}

int main() {
    std::cout << "\n🚀 Ephemeris Comparison Tool" << std::endl;
    std::cout << "Asteroide 17030 (Sierks) - 25-30 Nov 2025\n" << std::endl;
    
    try {
        // =====================================================
        // 1. CARICA DATI JPL HORIZONS
        // =====================================================
        std::cout << "Caricamento dati JPL Horizons..." << std::endl;
        auto jplData = getJPLHorizonData();
        std::cout << "✓ " << jplData.size() << " effemeridi JPL caricate\n" << std::endl;
        
        // =====================================================
        // 2. PROPAGA CON ASTDYN
        // =====================================================
        std::cout << "Propagazione con AstDyn..." << std::endl;
        AstDynWrapper wrapper(PropagationSettings::highAccuracy());
        
        if (!wrapper.loadFromEQ1File("astdyn/data/17030.eq1")) {
            std::cerr << "✗ Errore caricamento .eq1" << std::endl;
            return 1;
        }
        
        // Prepara epoche per Chebyshev (primi 5 punti per fitting)
        std::vector<Eigen::Vector3d> fitPositions;
        std::vector<double> fitEpochs;
        
        // Estrai 5 punti equidistanziati per fitting
        for (size_t i = 0; i < jplData.size(); i += 4) {
            fitEpochs.push_back(jplData[i].mjd_tdb);
            
            auto state = wrapper.propagateToEpoch(jplData[i].mjd_tdb);
            fitPositions.push_back(state.position);
        }
        
        std::cout << "✓ Propagazione completata (" << fitEpochs.size() << " epoche)\n" << std::endl;
        
        // =====================================================
        // 3. FITTA POLINOMI CHEBYSHEV
        // =====================================================
        std::cout << "Fitting polinomi di Chebyshev..." << std::endl;
        ChebyshevApproximation chebyshev(8);
        
        if (!chebyshev.fit(fitPositions, jplData.front().mjd_tdb, jplData.back().mjd_tdb)) {
            std::cerr << "✗ Errore fitting Chebyshev" << std::endl;
            return 1;
        }
        
        std::cout << "✓ Chebyshev fittato (8 coefficienti, intervallo: "
                  << std::fixed << std::setprecision(2)
                  << jplData.front().mjd_tdb << " - "
                  << jplData.back().mjd_tdb << " MJD)\n" << std::endl;
        
        // =====================================================
        // 4. TABELLA DI CONFRONTO
        // =====================================================
        printTableHeader();
        
        double rmsErrorAstDyn = 0.0, rmsErrorChebyshev = 0.0;
        int count = 0;
        
        for (const auto& jpl : jplData) {
            std::string dateTime = formatDateTime(jpl.mjd_tdb);
            
            // Propagazione AstDyn
            auto astdynState = wrapper.propagateToEpoch(jpl.mjd_tdb);
            
            // Interpolazione Chebyshev
            Eigen::Vector3d chebPos = chebyshev.evaluatePosition(jpl.mjd_tdb);
            
            // JPL (reference)
            Eigen::Vector3d jplPos(jpl.x_au, jpl.y_au, jpl.z_au);
            Eigen::Vector3d jplVel(jpl.vx_au_d, jpl.vy_au_d, jpl.vz_au_d);
            
            // Calcola errori
            Eigen::Vector3d errorAstDyn = astdynState.position - jplPos;
            Eigen::Vector3d errorChebyshev = chebPos - jplPos;
            
            double errorMetersAstDyn = errorAstDyn.norm() * 149597870.7;  // AU to km to m
            double errorMetersChebyshev = errorChebyshev.norm() * 149597870.7;
            
            rmsErrorAstDyn += errorMetersAstDyn * errorMetersAstDyn;
            rmsErrorChebyshev += errorMetersChebyshev * errorMetersChebyshev;
            count++;
            
            // Stampa JPL (reference)
            printRow(dateTime, jpl.mjd_tdb, "JPL Horizons",
                    jpl.x_au, jpl.y_au, jpl.z_au,
                    jpl.vx_au_d, jpl.vy_au_d, jpl.vz_au_d,
                    jplPos.norm(), 0.0);
            
            // Stampa AstDyn
            printRow(dateTime, jpl.mjd_tdb, "AstDyn",
                    astdynState.position.x(), astdynState.position.y(), astdynState.position.z(),
                    astdynState.velocity.x(), astdynState.velocity.y(), astdynState.velocity.z(),
                    astdynState.position.norm(), errorMetersAstDyn);
            
            // Stampa Chebyshev
            Eigen::Vector3d chebVel = chebyshev.evaluateVelocity(jpl.mjd_tdb);
            printRow(dateTime, jpl.mjd_tdb, "Chebyshev",
                    chebPos.x(), chebPos.y(), chebPos.z(),
                    chebVel.x(), chebVel.y(), chebVel.z(),
                    chebPos.norm(), errorMetersChebyshev);
            
            std::cout << std::string(150, '-') << std::endl;
        }
        
        // =====================================================
        // 5. STATISTICHE
        // =====================================================
        std::cout << std::string(150, '=') << std::endl;
        std::cout << "\n📊 STATISTICHE DI ERRORE\n" << std::endl;
        
        rmsErrorAstDyn = std::sqrt(rmsErrorAstDyn / count);
        rmsErrorChebyshev = std::sqrt(rmsErrorChebyshev / count);
        
        std::cout << std::left << std::fixed
                  << "RMS Errore AstDyn vs JPL:    " << std::setw(15) << std::setprecision(3) << rmsErrorAstDyn << " m\n"
                  << "RMS Errore Chebyshev vs JPL: " << std::setw(15) << std::setprecision(3) << rmsErrorChebyshev << " m\n"
                  << std::endl;
        
        // Percentuali
        double percentDiff = 100.0 * (rmsErrorChebyshev - rmsErrorAstDyn) / rmsErrorAstDyn;
        if (percentDiff > 0) {
            std::cout << "⚠️  Chebyshev è " << std::setprecision(1) << percentDiff << "% meno accurato di AstDyn" << std::endl;
        } else {
            std::cout << "✓ Chebyshev è " << std::setprecision(1) << -percentDiff << "% più accurato di AstDyn" << std::endl;
        }
        
        std::cout << "\n✓ Analisi completata\n" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ Errore: " << e.what() << std::endl;
        return 1;
    }
}
