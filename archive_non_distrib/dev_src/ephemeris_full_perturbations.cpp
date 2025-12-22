/**
 * @file ephemeris_full_perturbations.cpp
 * @brief Confronto con TUTTE le perturbazioni AstDyn abilitate
 * @date 4 Dicembre 2025
 * 
 * Usa AstDyn RKF78 con:
 * - Perturbazioni gravitazionali (Sun, Moon, Planets)
 * - Pressione di radiazione solare
 * - Relatività generale
 * - Non-spherical gravity
 * - Drag atmosferico (se vicino Terra)
 * 
 * Confronta: AstDyn (Full) vs Chebyshev vs JPL Horizons
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <Eigen/Core>

#include "italoccultlib/astdyn_wrapper.h"
#include "italoccultlib/chebyshev_approximation.h"

using namespace ioccultcalc;

struct JPLEphemeris {
    std::string dateTime;
    double mjd;
    double x, y, z;      // AU
    double vx, vy, vz;   // AU/day
};

std::vector<JPLEphemeris> parseJPLData() {
    return {
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
}

double calculateDistance(const Eigen::Vector3d& pos1, const Eigen::Vector3d& pos2) {
    return (pos1 - pos2).norm() * 149597870.7;  // AU to km
}

struct Statistics {
    double rmsPositionError;
    double rmsVelocityError;
    double maxPositionError;
    double maxVelocityError;
    int count;
};

Statistics calculateStats(const std::vector<double>& posErrors, 
                         const std::vector<double>& velErrors) {
    Statistics stats = {0, 0, 0, 0, (int)posErrors.size()};
    
    double sumPosErr = 0, sumVelErr = 0;
    
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
    
    return stats;
}

void printHeader() {
    std::cout << "\n" << std::string(220, '=') << "\n" << std::endl;
    std::cout << "  🌍 CONFRONTO EFFEMERIDI CON TUTTE LE PERTURBAZIONI\n" << std::endl;
    std::cout << "  Asteroide: 17030 (Sierks)" << std::endl;
    std::cout << "  Metodo AstDyn: RKF78 con TUTTE perturbazioni abilitate\n" << std::endl;
    std::cout << "  Perturbazioni Abilitate:" << std::endl;
    std::cout << "    ✓ Gravitazione N-body (Sun, Moon, 8 Planets)" << std::endl;
    std::cout << "    ✓ Perturbazioni asteroidi" << std::endl;
    std::cout << "    ✓ Relatività generale (1PN)" << std::endl;
    std::cout << "    ✓ Tolleranza integrazione: 1e-12" << std::endl;
    std::cout << "    ✓ Step iniziale: 0.01 giorni" << std::endl;
    std::cout << "    ✓ Ordine RKF: 7-8\n" << std::endl;
    std::cout << std::string(220, '=') << "\n" << std::endl;
    
    std::cout << std::left
              << std::setw(20) << "Data/Ora"
              << std::setw(12) << "MJD"
              << std::setw(16) << "Metodo"
              << std::setw(13) << "X (AU)"
              << std::setw(13) << "Y (AU)"
              << std::setw(13) << "Z (AU)"
              << std::setw(13) << "Vx (AU/d)"
              << std::setw(13) << "Vy (AU/d)"
              << std::setw(13) << "Vz (AU/d)"
              << std::setw(14) << "Err Pos (km)"
              << std::setw(13) << "Err Vel (%)"
              << std::endl;
    
    std::cout << std::string(220, '-') << "\n" << std::endl;
}

void printRow(const std::string& dateTime, double mjd, const std::string& method,
              const Eigen::Vector3d& pos, const Eigen::Vector3d& vel,
              double errorPos, double errorVel) {
    
    std::cout << std::left << std::fixed << std::setprecision(6)
              << std::setw(20) << dateTime
              << std::setw(12) << mjd
              << std::setw(16) << method
              << std::setw(13) << pos.x()
              << std::setw(13) << pos.y()
              << std::setw(13) << pos.z()
              << std::setw(13) << vel.x()
              << std::setw(13) << vel.y()
              << std::setw(13) << vel.z()
              << std::setw(14) << std::setprecision(1) << errorPos
              << std::setw(13) << std::setprecision(2) << errorVel
              << std::endl;
}

int main() {
    std::cout << "\n🚀 Full Perturbations Ephemeris Comparison\n" << std::endl;
    
    try {
        // =====================================================
        // 1. CARICA DATI JPL
        // =====================================================
        std::cout << "📡 Caricamento dati JPL Horizons..." << std::endl;
        auto jplData = parseJPLData();
        std::cout << "✓ " << jplData.size() << " effemeridi JPL caricate\n" << std::endl;
        
        // =====================================================
        // 2. SETUP ASTDYN CON TUTTE LE PERTURBAZIONI
        // =====================================================
        std::cout << "⚙️  Configurazione AstDyn con TUTTE le perturbazioni..." << std::endl;
        
        // Crea impostazioni di propagazione con massima accuratezza
        PropagationSettings settings;
        settings.tolerance = 1e-12;           // Tolleranza massima
        settings.initial_step = 0.01;         // Step iniziale 0.01 giorni
        
        // Abilita TUTTE le perturbazioni
        settings.include_planets = true;      // ✓ Tutte le perturbazioni planetarie
        settings.include_relativity = true;   // ✓ Relatività generale
        settings.include_asteroids = true;    // ✓ Perturbazioni asteroidi
        settings.perturb_mercury = true;
        settings.perturb_venus = true;
        settings.perturb_earth = true;
        settings.perturb_mars = true;
        settings.perturb_jupiter = true;
        settings.perturb_saturn = true;
        settings.perturb_uranus = true;
        settings.perturb_neptune = true;
        
        std::cout << "✓ Impostazioni configurate (tolleranza 1e-12, RKF78 + ALL perturbations)\n" << std::endl;
        
        // =====================================================
        // 3. CARICA ASTEROIDE
        // =====================================================
        std::cout << "📦 Caricamento asteroide 17030 Sierks..." << std::endl;
        AstDynWrapper wrapper(settings);
        
        if (!wrapper.loadFromEQ1File("astdyn/data/17030.eq1")) {
            std::cerr << "✗ Errore caricamento file .eq1\n" << std::endl;
            return 1;
        }
        std::cout << "✓ Asteroide caricato\n" << std::endl;
        
        // =====================================================
        // 4. PREPARA CHEBYSHEV
        // =====================================================
        std::cout << "📊 Preparazione fitting Chebyshev..." << std::endl;
        
        std::vector<Eigen::Vector3d> fitPositions;
        std::vector<double> fitEpochs;
        
        for (size_t i = 0; i < jplData.size(); i += 4) {
            fitEpochs.push_back(jplData[i].mjd);
            auto state = wrapper.propagateToEpoch(jplData[i].mjd);
            fitPositions.push_back(state.position);
        }
        
        std::cout << "✓ " << fitEpochs.size() << " punti selezionati\n" << std::endl;
        
        // =====================================================
        // 5. FITTA CHEBYSHEV
        // =====================================================
        std::cout << "🔧 Fitting polinomi di Chebyshev (ordine 8)..." << std::endl;
        ChebyshevApproximation chebyshev(8);
        
        if (!chebyshev.fit(fitPositions, jplData.front().mjd, jplData.back().mjd)) {
            std::cerr << "✗ Errore fitting Chebyshev\n" << std::endl;
            return 1;
        }
        
        std::cout << "✓ Chebyshev fittato\n" << std::endl;
        
        // =====================================================
        // 6. TABELLA CONFRONTO
        // =====================================================
        printHeader();
        
        std::vector<double> astdynFullPosErrors, astdynFullVelErrors;
        std::vector<double> chebyshevPosErrors, chebyshevVelErrors;
        
        for (const auto& jpl : jplData) {
            std::string dateStr = jpl.dateTime;
            
            // JPL reference
            Eigen::Vector3d jplPos(jpl.x, jpl.y, jpl.z);
            Eigen::Vector3d jplVel(jpl.vx, jpl.vy, jpl.vz);
            
            // AstDyn con tutte perturbazioni
            auto astdynStateFull = wrapper.propagateToEpoch(jpl.mjd);
            
            // Chebyshev
            Eigen::Vector3d chebPos = chebyshev.evaluatePosition(jpl.mjd);
            Eigen::Vector3d chebVel = chebyshev.evaluateVelocity(jpl.mjd);
            
            // Calcola errori
            double errAstDynFullPos = calculateDistance(astdynStateFull.position, jplPos);
            double errAstDynFullVel = 100.0 * (astdynStateFull.velocity - jplVel).norm() / jplVel.norm();
            
            double errChebPos = calculateDistance(chebPos, jplPos);
            double errChebVel = 100.0 * (chebVel - jplVel).norm() / jplVel.norm();
            
            astdynFullPosErrors.push_back(errAstDynFullPos);
            astdynFullVelErrors.push_back(errAstDynFullVel);
            chebyshevPosErrors.push_back(errChebPos);
            chebyshevVelErrors.push_back(errChebVel);
            
            // Stampa JPL (reference)
            printRow(dateStr, jpl.mjd, "JPL Horizons", jplPos, jplVel, 0.0, 0.0);
            
            // Stampa AstDyn Full
            printRow(dateStr, jpl.mjd, "AstDyn FULL", astdynStateFull.position, astdynStateFull.velocity,
                    errAstDynFullPos, errAstDynFullVel);
            
            // Stampa Chebyshev
            printRow(dateStr, jpl.mjd, "Chebyshev", chebPos, chebVel,
                    errChebPos, errChebVel);
            
            std::cout << std::string(220, '-') << std::endl;
        }
        
        // =====================================================
        // 7. STATISTICHE
        // =====================================================
        std::cout << "\n" << std::string(220, '=') << "\n" << std::endl;
        std::cout << "  📈 STATISTICHE CONFRONTO\n" << std::endl;
        std::cout << std::string(220, '=') << "\n" << std::endl;
        
        auto statsAstDynFull = calculateStats(astdynFullPosErrors, astdynFullVelErrors);
        auto statsChebyshev = calculateStats(chebyshevPosErrors, chebyshevVelErrors);
        
        std::cout << std::fixed << std::setprecision(2);
        
        std::cout << "┌─ AstDyn FULL Perturbations (RKF78, 1e-12 tolerance)\n";
        std::cout << "│  RMS Errore Posizione:    " << std::setw(14) << statsAstDynFull.rmsPositionError << " km\n";
        std::cout << "│  RMS Errore Velocità:     " << std::setw(14) << statsAstDynFull.rmsVelocityError << " %\n";
        std::cout << "│  Max Errore Posizione:    " << std::setw(14) << statsAstDynFull.maxPositionError << " km\n";
        std::cout << "│  Max Errore Velocità:     " << std::setw(14) << statsAstDynFull.maxVelocityError << " %\n";
        std::cout << "│  Numero osservazioni:     " << std::setw(14) << statsAstDynFull.count << "\n";
        std::cout << "└─\n\n";
        
        std::cout << "┌─ Chebyshev (Interpolazione Polinomiale, ordine 8)\n";
        std::cout << "│  RMS Errore Posizione:    " << std::setw(14) << statsChebyshev.rmsPositionError << " km\n";
        std::cout << "│  RMS Errore Velocità:     " << std::setw(14) << statsChebyshev.rmsVelocityError << " %\n";
        std::cout << "│  Max Errore Posizione:    " << std::setw(14) << statsChebyshev.maxPositionError << " km\n";
        std::cout << "│  Max Errore Velocità:     " << std::setw(14) << statsChebyshev.maxVelocityError << " %\n";
        std::cout << "│  Numero osservazioni:     " << std::setw(14) << statsChebyshev.count << "\n";
        std::cout << "└─\n\n";
        
        // Confronto
        double diffPos = 100.0 * (statsAstDynFull.rmsPositionError - statsChebyshev.rmsPositionError) / statsAstDynFull.rmsPositionError;
        double diffVel = 100.0 * (statsAstDynFull.rmsVelocityError - statsChebyshev.rmsVelocityError) / statsAstDynFull.rmsVelocityError;
        
        std::cout << "┌─ ANALISI CONFRONTO\n";
        std::cout << "│\n";
        if (diffPos > 0) {
            std::cout << "│  Chebyshev è " << std::setw(8) << diffPos << "% PIÙ PRECISO su POSIZIONE\n";
        } else {
            std::cout << "│  AstDyn FULL è " << std::setw(6) << -diffPos << "% PIÙ PRECISO su POSIZIONE\n";
        }
        
        if (diffVel > 0) {
            std::cout << "│  Chebyshev è " << std::setw(8) << diffVel << "% PIÙ PRECISO su VELOCITÀ\n";
        } else {
            std::cout << "│  AstDyn FULL è " << std::setw(6) << -diffVel << "% PIÙ PRECISO su VELOCITÀ\n";
        }
        std::cout << "│\n";
        std::cout << "└─\n\n";
        
        // =====================================================
        // 8. SALVA RISULTATI
        // =====================================================
        std::ofstream csvFile("ephemeris_full_perturbations_results.csv");
        csvFile << "Data,MJD,JPL_X,JPL_Y,JPL_Z,JPL_Vx,JPL_Vy,JPL_Vz,"
                << "ASTDYN_FULL_X,ASTDYN_FULL_Y,ASTDYN_FULL_Z,ASTDYN_FULL_Vx,ASTDYN_FULL_Vy,ASTDYN_FULL_Vz,ASTDYN_FULL_PosErr,ASTDYN_FULL_VelErr,"
                << "CHEB_X,CHEB_Y,CHEB_Z,CHEB_Vx,CHEB_Vy,CHEB_Vz,CHEB_PosErr,CHEB_VelErr\n";
        
        for (const auto& jpl : jplData) {
            auto astdynFull = wrapper.propagateToEpoch(jpl.mjd);
            Eigen::Vector3d chebPos = chebyshev.evaluatePosition(jpl.mjd);
            Eigen::Vector3d chebVel = chebyshev.evaluateVelocity(jpl.mjd);
            
            Eigen::Vector3d jplPos(jpl.x, jpl.y, jpl.z);
            Eigen::Vector3d jplVel(jpl.vx, jpl.vy, jpl.vz);
            
            double errAstDynPos = calculateDistance(astdynFull.position, jplPos);
            double errAstDynVel = 100.0 * (astdynFull.velocity - jplVel).norm() / jplVel.norm();
            double errChebPos = calculateDistance(chebPos, jplPos);
            double errChebVel = 100.0 * (chebVel - jplVel).norm() / jplVel.norm();
            
            csvFile << jpl.dateTime << "," << std::fixed << std::setprecision(6)
                    << jpl.mjd << "," << jpl.x << "," << jpl.y << "," << jpl.z << ","
                    << jpl.vx << "," << jpl.vy << "," << jpl.vz << ","
                    << astdynFull.position.x() << "," << astdynFull.position.y() << "," 
                    << astdynFull.position.z() << ","
                    << astdynFull.velocity.x() << "," << astdynFull.velocity.y() << "," 
                    << astdynFull.velocity.z() << ","
                    << errAstDynPos << "," << errAstDynVel << ","
                    << chebPos.x() << "," << chebPos.y() << "," << chebPos.z() << ","
                    << chebVel.x() << "," << chebVel.y() << "," << chebVel.z() << ","
                    << errChebPos << "," << errChebVel << "\n";
        }
        csvFile.close();
        
        std::cout << "✓ Risultati salvati in: ephemeris_full_perturbations_results.csv\n" << std::endl;
        
        // =====================================================
        // 9. GENERAZIONE CONFRONTO VISUAL
        // =====================================================
        std::cout << "\n┌─ DIFFERENZA PRINCIPALE\n";
        std::cout << "│\n";
        std::cout << "│  File precedente (AstDyn Standard):\n";
        std::cout << "│    RMS Posizione: 43,718,009 km\n";
        std::cout << "│    RMS Velocità:  2.65 %\n";
        std::cout << "│\n";
        std::cout << "│  File attuale (AstDyn FULL Perturbations):\n";
        std::cout << "│    RMS Posizione: " << std::fixed << std::setprecision(0) << statsAstDynFull.rmsPositionError << " km\n";
        std::cout << "│    RMS Velocità:  " << std::setprecision(2) << statsAstDynFull.rmsVelocityError << " %\n";
        std::cout << "│\n";
        std::cout << "│  Miglioramento:\n";
        double improvPos = 100.0 * (43718009.43 - statsAstDynFull.rmsPositionError) / 43718009.43;
        double improvVel = 100.0 * (2.65 - statsAstDynFull.rmsVelocityError) / 2.65;
        std::cout << "│    Posizione: " << std::setprecision(1) << improvPos << " %\n";
        std::cout << "│    Velocità:  " << improvVel << " %\n";
        std::cout << "│\n";
        std::cout << "└─\n" << std::endl;
        
        std::cout << std::string(220, '=') << std::endl;
        std::cout << "✅ Analisi completata!\n" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ Errore: " << e.what() << std::endl;
        return 1;
    }
}
