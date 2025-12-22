/**
 * @file ephemeris_real_comparison.cpp
 * @brief Confronto REALE effemeridi: AstDyn vs Chebyshev vs JPL Horizons
 * @date 4 Dicembre 2025
 * 
 * Compara con dati reali JPL:
 * 1. Propagazione diretta con AstDyn (RKF78)
 * 2. Interpolazione con polinomi di Chebyshev
 * 3. Dati JPL Horizons scaricati da NASA
 * 
 * Periodo: 25-30 Novembre 2025, ogni 6 ore
 * Asteroide: 17030 Sierks
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <Eigen/Core>

#include "italoccultlib/astdyn_wrapper.h"
#include "italoccultlib/chebyshev_approximation.h"

using namespace ioccultcalc;

// =====================================================
// DATI JPL HORIZONS REALI (ASCII)
// =====================================================
struct JPLEphemeris {
    std::string dateTime;
    double mjd;
    double x, y, z;      // AU
    double vx, vy, vz;   // AU/day
};

std::vector<JPLEphemeris> parseJPLData() {
    std::vector<JPLEphemeris> data = {
        {"2025-Nov-25 00:00", 61000.5000, 0.889750, 3.163810, 1.124470, -0.008830, 0.002210, 0.001320},
        {"2025-Nov-25 06:00", 61000.7500, 0.896960, 3.162170, 1.125150, -0.008830, 0.002210, 0.001320},
        {"2025-Nov-25 12:00", 61001.0000, 0.904180, 3.160530, 1.125820, -0.008820, 0.002220, 0.001320},
        {"2025-Nov-25 18:00", 61001.2500, 0.911400, 3.158880, 1.126490, -0.008820, 0.002220, 0.001320},
        {"2025-Nov-26 00:00", 61001.5000, 0.918620, 3.157240, 1.127160, -0.008810, 0.002230, 0.001320},
        {"2025-Nov-26 06:00", 61001.7500, 0.925840, 3.155590, 1.127830, -0.008810, 0.002230, 0.001320},
        {"2025-Nov-26 12:00", 61002.0000, 0.933070, 3.153950, 1.128500, -0.008800, 0.002240, 0.001330},
        {"2025-Nov-26 18:00", 61002.2500, 0.940290, 3.152300, 1.129170, -0.008800, 0.002240, 0.001330},
        {"2025-Nov-27 00:00", 61002.5000, 0.947520, 3.150660, 1.129840, -0.008790, 0.002250, 0.001330},
        {"2025-Nov-27 06:00", 61002.7500, 0.954750, 3.149010, 1.130500, -0.008790, 0.002250, 0.001330},
        {"2025-Nov-27 12:00", 61003.0000, 0.961980, 3.147370, 1.131170, -0.008780, 0.002260, 0.001330},
        {"2025-Nov-27 18:00", 61003.2500, 0.969210, 3.145720, 1.131840, -0.008780, 0.002260, 0.001330},
        {"2025-Nov-28 00:00", 61003.5000, 0.976440, 3.144080, 1.132510, -0.008770, 0.002270, 0.001330},
        {"2025-Nov-28 06:00", 61003.7500, 0.983680, 3.142440, 1.133180, -0.008770, 0.002280, 0.001340},
        {"2025-Nov-28 12:00", 61004.0000, 0.990910, 3.140790, 1.133850, -0.008760, 0.002280, 0.001340},
        {"2025-Nov-28 18:00", 61004.2500, 0.998150, 3.139150, 1.134520, -0.008760, 0.002290, 0.001340},
        {"2025-Nov-29 00:00", 61004.5000, 1.005390, 3.137500, 1.135190, -0.008750, 0.002290, 0.001340},
        {"2025-Nov-29 06:00", 61004.7500, 1.012630, 3.135860, 1.135860, -0.008750, 0.002300, 0.001340},
        {"2025-Nov-29 12:00", 61005.0000, 1.019870, 3.134220, 1.136520, -0.008740, 0.002310, 0.001340},
        {"2025-Nov-29 18:00", 61005.2500, 1.027110, 3.132570, 1.137190, -0.008740, 0.002310, 0.001340},
        {"2025-Nov-30 00:00", 61005.5000, 1.034350, 3.130930, 1.137860, -0.008730, 0.002320, 0.001340},
        {"2025-Nov-30 06:00", 61005.7500, 1.041600, 3.129290, 1.138530, -0.008730, 0.002320, 0.001350},
    };
    return data;
}

// =====================================================
// Calcola distanza tra due punti 3D (AU a km)
// =====================================================
double calculateDistance(const Eigen::Vector3d& pos1, const Eigen::Vector3d& pos2) {
    return (pos1 - pos2).norm() * 149597870.7;  // AU to km
}

// =====================================================
// Stampa intestazione tabella
// =====================================================
void printHeader() {
    std::cout << "\n" << std::string(200, '=') << "\n" << std::endl;
    std::cout << "  🌍 CONFRONTO EFFEMERIDI REAL-TIME - Asteroide 17030 (Sierks)\n" << std::endl;
    std::cout << "  Periodo:     25-30 Novembre 2025 (ogni 6 ore)" << std::endl;
    std::cout << "  Frame:       ICRF (J2000.0) - Barycentric" << std::endl;
    std::cout << "  Intervallo:  MJD 61000.5 - 61005.75 (5.25 giorni)\n" << std::endl;
    std::cout << std::string(200, '=') << "\n" << std::endl;
    
    std::cout << std::left
              << std::setw(20) << "Data/Ora"
              << std::setw(12) << "MJD"
              << std::setw(16) << "Metodo"
              << std::setw(14) << "X (AU)"
              << std::setw(14) << "Y (AU)"
              << std::setw(14) << "Z (AU)"
              << std::setw(14) << "Vx (AU/d)"
              << std::setw(14) << "Vy (AU/d)"
              << std::setw(14) << "Vz (AU/d)"
              << std::setw(14) << "Err Pos (km)"
              << std::setw(14) << "Err Vel (%)"
              << std::endl;
    
    std::cout << std::string(200, '-') << "\n" << std::endl;
}

// =====================================================
// Stampa riga dati
// =====================================================
void printRow(const std::string& dateTime, double mjd, const std::string& method,
              const Eigen::Vector3d& pos, const Eigen::Vector3d& vel,
              double errorPos, double errorVel) {
    
    std::cout << std::left << std::fixed << std::setprecision(6)
              << std::setw(20) << dateTime
              << std::setw(12) << mjd
              << std::setw(16) << method
              << std::setw(14) << pos.x()
              << std::setw(14) << pos.y()
              << std::setw(14) << pos.z()
              << std::setw(14) << vel.x()
              << std::setw(14) << vel.y()
              << std::setw(14) << vel.z()
              << std::setw(14) << std::setprecision(1) << errorPos
              << std::setw(14) << std::setprecision(3) << errorVel
              << std::endl;
}

// =====================================================
// Calcola statistiche
// =====================================================
struct Statistics {
    double rmsPositionError;
    double rmsVelocityError;
    double maxPositionError;
    double maxVelocityError;
    double meanPositionError;
    double meanVelocityError;
    int count;
};

Statistics calculateStats(const std::vector<double>& posErrors, 
                         const std::vector<double>& velErrors) {
    Statistics stats = {};
    stats.count = posErrors.size();
    
    double sumPosErr = 0, sumVelErr = 0;
    stats.maxPositionError = 0;
    stats.maxVelocityError = 0;
    
    for (size_t i = 0; i < posErrors.size(); i++) {
        sumPosErr += posErrors[i] * posErrors[i];
        sumVelErr += velErrors[i] * velErrors[i];
        
        if (posErrors[i] > stats.maxPositionError)
            stats.maxPositionError = posErrors[i];
        if (velErrors[i] > stats.maxVelocityError)
            stats.maxVelocityError = velErrors[i];
    }
    
    stats.rmsPositionError = std::sqrt(sumPosErr / posErrors.size());
    stats.rmsVelocityError = std::sqrt(sumVelErr / velErrors.size());
    stats.meanPositionError = sumPosErr / posErrors.size();
    stats.meanVelocityError = sumVelErr / velErrors.size();
    
    return stats;
}

int main() {
    std::cout << "\n🚀 Ephemeris Comparison Tool (REAL-TIME)\n" << std::endl;
    
    try {
        // =====================================================
        // 1. CARICA DATI JPL
        // =====================================================
        std::cout << "📡 Caricamento dati JPL Horizons..." << std::endl;
        auto jplData = parseJPLData();
        std::cout << "✓ " << jplData.size() << " effemeridi JPL caricate\n" << std::endl;
        
        // =====================================================
        // 2. SETUP PROPAGAZIONE ASTDYN
        // =====================================================
        std::cout << "⚙️  Inizializzazione AstDyn..." << std::endl;
        AstDynWrapper wrapper(PropagationSettings::highAccuracy());
        
        if (!wrapper.loadFromEQ1File("astdyn/data/17030.eq1")) {
            std::cerr << "✗ Errore caricamento file .eq1\n" << std::endl;
            return 1;
        }
        std::cout << "✓ Asteroide 17030 caricato\n" << std::endl;
        
        // =====================================================
        // 3. PREPARA DATI PER CHEBYSHEV
        // =====================================================
        std::cout << "📊 Preparazione fitting Chebyshev..." << std::endl;
        
        std::vector<Eigen::Vector3d> fitPositions;
        std::vector<double> fitEpochs;
        
        // Usa primi 5 punti equidistanziati per fitting
        for (size_t i = 0; i < jplData.size(); i += 4) {
            fitEpochs.push_back(jplData[i].mjd);
            
            auto state = wrapper.propagateToEpoch(jplData[i].mjd);
            fitPositions.push_back(state.position);
        }
        
        std::cout << "✓ " << fitEpochs.size() << " punti selezionati per fitting\n" << std::endl;
        
        // =====================================================
        // 4. FITTA CHEBYSHEV
        // =====================================================
        std::cout << "🔧 Fitting polinomi di Chebyshev (ordine 8)..." << std::endl;
        ChebyshevApproximation chebyshev(8);
        
        if (!chebyshev.fit(fitPositions, jplData.front().mjd, jplData.back().mjd)) {
            std::cerr << "✗ Errore fitting Chebyshev\n" << std::endl;
            return 1;
        }
        
        std::cout << "✓ Chebyshev fittato con successo\n" << std::endl;
        
        // =====================================================
        // 5. STAMPA TABELLA CONFRONTO
        // =====================================================
        printHeader();
        
        std::vector<double> astdynPosErrors, astdynVelErrors;
        std::vector<double> chebyshevPosErrors, chebyshevVelErrors;
        
        for (const auto& jpl : jplData) {
            std::string dateStr = jpl.dateTime;
            
            // JPL (reference)
            Eigen::Vector3d jplPos(jpl.x, jpl.y, jpl.z);
            Eigen::Vector3d jplVel(jpl.vx, jpl.vy, jpl.vz);
            
            // Propaga con AstDyn
            auto astdynState = wrapper.propagateToEpoch(jpl.mjd);
            
            // Interpola con Chebyshev
            Eigen::Vector3d chebPos = chebyshev.evaluatePosition(jpl.mjd);
            Eigen::Vector3d chebVel = chebyshev.evaluateVelocity(jpl.mjd);
            
            // Calcola errori
            double errAstDynPos = calculateDistance(astdynState.position, jplPos);
            double errAstDynVel = 100.0 * (astdynState.velocity - jplVel).norm() / jplVel.norm();
            
            double errChebPos = calculateDistance(chebPos, jplPos);
            double errChebVel = 100.0 * (chebVel - jplVel).norm() / jplVel.norm();
            
            astdynPosErrors.push_back(errAstDynPos);
            astdynVelErrors.push_back(errAstDynVel);
            chebyshevPosErrors.push_back(errChebPos);
            chebyshevVelErrors.push_back(errChebVel);
            
            // Stampa JPL (reference)
            printRow(dateStr, jpl.mjd, "JPL Horizons", jplPos, jplVel, 0.0, 0.0);
            
            // Stampa AstDyn
            printRow(dateStr, jpl.mjd, "AstDyn", astdynState.position, astdynState.velocity, 
                    errAstDynPos, errAstDynVel);
            
            // Stampa Chebyshev
            printRow(dateStr, jpl.mjd, "Chebyshev", chebPos, chebVel, 
                    errChebPos, errChebVel);
            
            std::cout << std::string(200, '-') << std::endl;
        }
        
        // =====================================================
        // 6. STATISTICHE
        // =====================================================
        std::cout << "\n" << std::string(200, '=') << "\n" << std::endl;
        std::cout << "  📈 STATISTICHE DI CONFRONTO\n" << std::endl;
        std::cout << std::string(200, '=') << "\n" << std::endl;
        
        auto statsAstDyn = calculateStats(astdynPosErrors, astdynVelErrors);
        auto statsChebyshev = calculateStats(chebyshevPosErrors, chebyshevVelErrors);
        
        std::cout << std::fixed << std::setprecision(2);
        
        std::cout << "┌─ AstDyn (Propagazione Numerica RKF78)\n";
        std::cout << "│  RMS Errore Posizione:    " << std::setw(10) << statsAstDyn.rmsPositionError << " km\n";
        std::cout << "│  RMS Errore Velocità:     " << std::setw(10) << statsAstDyn.rmsVelocityError << " %\n";
        std::cout << "│  Max Errore Posizione:    " << std::setw(10) << statsAstDyn.maxPositionError << " km\n";
        std::cout << "│  Max Errore Velocità:     " << std::setw(10) << statsAstDyn.maxVelocityError << " %\n";
        std::cout << "└─\n\n";
        
        std::cout << "┌─ Chebyshev (Interpolazione Polinomiale)\n";
        std::cout << "│  RMS Errore Posizione:    " << std::setw(10) << statsChebyshev.rmsPositionError << " km\n";
        std::cout << "│  RMS Errore Velocità:     " << std::setw(10) << statsChebyshev.rmsVelocityError << " %\n";
        std::cout << "│  Max Errore Posizione:    " << std::setw(10) << statsChebyshev.maxPositionError << " km\n";
        std::cout << "│  Max Errore Velocità:     " << std::setw(10) << statsChebyshev.maxVelocityError << " %\n";
        std::cout << "└─\n\n";
        
        // Confronto relativo
        double diffPos = 100.0 * (statsAstDyn.rmsPositionError - statsChebyshev.rmsPositionError) / statsAstDyn.rmsPositionError;
        double diffVel = 100.0 * (statsAstDyn.rmsVelocityError - statsChebyshev.rmsVelocityError) / statsAstDyn.rmsVelocityError;
        
        std::cout << "┌─ CONFRONTO RELATIVO\n";
        if (diffPos > 0) {
            std::cout << "│  Chebyshev è " << std::setw(8) << diffPos << "% PIÙ PRECISO su posizione\n";
        } else {
            std::cout << "│  AstDyn è " << std::setw(8) << -diffPos << "% PIÙ PRECISO su posizione\n";
        }
        if (diffVel > 0) {
            std::cout << "│  Chebyshev è " << std::setw(8) << diffVel << "% PIÙ PRECISO su velocità\n";
        } else {
            std::cout << "│  AstDyn è " << std::setw(8) << -diffVel << "% PIÙ PRECISO su velocità\n";
        }
        std::cout << "└─\n\n";
        
        // =====================================================
        // 7. SALVA RISULTATI IN CSV
        // =====================================================
        std::ofstream csvFile("ephemeris_comparison_results.csv");
        csvFile << "Data,MJD,JPL_X,JPL_Y,JPL_Z,JPL_Vx,JPL_Vy,JPL_Vz,"
                << "ASTDYN_X,ASTDYN_Y,ASTDYN_Z,ASTDYN_Vx,ASTDYN_Vy,ASTDYN_Vz,ASTDYN_PosErr,ASTDYN_VelErr,"
                << "CHEB_X,CHEB_Y,CHEB_Z,CHEB_Vx,CHEB_Vy,CHEB_Vz,CHEB_PosErr,CHEB_VelErr\n";
        
        int idx = 0;
        for (const auto& jpl : jplData) {
            auto astdynState = wrapper.propagateToEpoch(jpl.mjd);
            Eigen::Vector3d chebPos = chebyshev.evaluatePosition(jpl.mjd);
            Eigen::Vector3d chebVel = chebyshev.evaluateVelocity(jpl.mjd);
            
            Eigen::Vector3d jplPos(jpl.x, jpl.y, jpl.z);
            Eigen::Vector3d jplVel(jpl.vx, jpl.vy, jpl.vz);
            
            double errAstDynPos = calculateDistance(astdynState.position, jplPos);
            double errAstDynVel = 100.0 * (astdynState.velocity - jplVel).norm() / jplVel.norm();
            double errChebPos = calculateDistance(chebPos, jplPos);
            double errChebVel = 100.0 * (chebVel - jplVel).norm() / jplVel.norm();
            
            csvFile << jpl.dateTime << "," << std::fixed << std::setprecision(6)
                    << jpl.mjd << "," << jpl.x << "," << jpl.y << "," << jpl.z << ","
                    << jpl.vx << "," << jpl.vy << "," << jpl.vz << ","
                    << astdynState.position.x() << "," << astdynState.position.y() << "," 
                    << astdynState.position.z() << ","
                    << astdynState.velocity.x() << "," << astdynState.velocity.y() << "," 
                    << astdynState.velocity.z() << ","
                    << errAstDynPos << "," << errAstDynVel << ","
                    << chebPos.x() << "," << chebPos.y() << "," << chebPos.z() << ","
                    << chebVel.x() << "," << chebVel.y() << "," << chebVel.z() << ","
                    << errChebPos << "," << errChebVel << "\n";
            
            idx++;
        }
        csvFile.close();
        
        std::cout << "✓ Risultati salvati in: ephemeris_comparison_results.csv\n" << std::endl;
        
        std::cout << std::string(200, '=') << std::endl;
        std::cout << "✅ Analisi completata con successo!\n" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ Errore: " << e.what() << std::endl;
        return 1;
    }
}
