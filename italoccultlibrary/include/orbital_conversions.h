/**
 * @file orbital_conversions.h
 * @brief Conversioni tra sistemi di elementi orbitali
 * @author IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 * 
 * Questo modulo fornisce conversioni tra:
 * - Elementi equinoziali (AstDyS) → Elementi kepleriani
 * - Elementi kepleriani → Stato cartesiano (eclittico)
 * - Stato eclittico → Stato ICRF (equatoriale J2000)
 * 
 * Include anche:
 * - Risoluzione equazione di Keplero (M = E - e·sin(E))
 * - Validazione frame di riferimento
 */

#ifndef IOCCULTCALC_ORBITAL_CONVERSIONS_H
#define IOCCULTCALC_ORBITAL_CONVERSIONS_H

#include "eq1_parser.h"
#include <Eigen/Dense>
#include <cmath>

namespace ioccultcalc {

/**
 * @struct KeplerianElements
 * @brief Elementi orbitali kepleriani classici
 * 
 * Sistema kepleriano standard:
 * - a: semiasse maggiore [AU]
 * - e: eccentricità [adimensionale]
 * - i: inclinazione [radianti]
 * - Ω: longitudine nodo ascendente [radianti]
 * - ω: argomento del perielio [radianti]
 * - M: anomalia media [radianti]
 * 
 * Frame: ECLM J2000 (eclittica media J2000)
 */
struct KeplerianElements {
    double a;          ///< Semiasse maggiore [AU]
    double e;          ///< Eccentricità
    double i;          ///< Inclinazione [rad]
    double Omega;      ///< Longitudine nodo ascendente [rad]
    double omega;      ///< Argomento perielio [rad]
    double M;          ///< Anomalia media [rad]
    double epoch_jd;   ///< Epoca [Julian Date, TDT]
    std::string name;  ///< Nome oggetto
    
    /// Costruttore di default
    KeplerianElements()
        : a(0.0), e(0.0), i(0.0), Omega(0.0), omega(0.0), M(0.0),
          epoch_jd(0.0), name("") {}
    
    /// Validazione
    bool isValid() const {
        return a > 0.0 && e >= 0.0 && e < 1.0 && epoch_jd > 0.0;
    }
    
    /// Converti angoli da radianti a gradi
    void toDegreesInPlace() {
        constexpr double RAD_TO_DEG = 180.0 / M_PI;
        i *= RAD_TO_DEG;
        Omega *= RAD_TO_DEG;
        omega *= RAD_TO_DEG;
        M *= RAD_TO_DEG;
    }
};

/**
 * @struct CartesianState
 * @brief Stato orbitale cartesiano (posizione + velocità)
 * 
 * Rappresenta lo stato di un corpo celeste in coordinate cartesiane:
 * - position: vettore posizione [AU]
 * - velocity: vettore velocità [AU/day]
 * - epoch_jd: epoca dello stato [Julian Date]
 */
struct CartesianState {
    Eigen::Vector3d position;  ///< Posizione [AU]
    Eigen::Vector3d velocity;  ///< Velocità [AU/day]
    double epoch_jd;           ///< Epoca [JD]
    
    /// Costruttore di default
    CartesianState() 
        : position(Eigen::Vector3d::Zero()),
          velocity(Eigen::Vector3d::Zero()),
          epoch_jd(0.0) {}
    
    /// Costruttore con valori
    CartesianState(const Eigen::Vector3d& pos, 
                   const Eigen::Vector3d& vel,
                   double jd)
        : position(pos), velocity(vel), epoch_jd(jd) {}
    
    /// Energia orbitale specifica
    double getSpecificEnergy(double mu = 0.0002959122082855911) const {
        double r = position.norm();
        double v2 = velocity.squaredNorm();
        return 0.5 * v2 - mu / r;
    }
    
    /// Momento angolare specifico
    Eigen::Vector3d getAngularMomentum() const {
        return position.cross(velocity);
    }
};

/**
 * @class OrbitalConversions
 * @brief Classe statica per conversioni tra sistemi orbitali
 * 
 * Esempio d'uso:
 * @code
 *   // Parse elementi da .eq1
 *   auto eq = EQ1Parser::parseFile("17030.eq1");
 *   
 *   // Converti equinoziali → kepleriani
 *   auto kep = OrbitalConversions::equinoctialToKeplerian(eq);
 *   
 *   // Converti kepleriani → cartesiano (eclittico)
 *   auto state_ecl = OrbitalConversions::keplerianToCartesian(kep);
 *   
 *   // Converti eclittico → ICRF
 *   auto state_icrf = OrbitalConversions::eclipticToICRF(state_ecl);
 * @endcode
 */
class OrbitalConversions {
public:
    // ========================================================================
    // Costanti astronomiche
    // ========================================================================
    
    static constexpr double PI = 3.141592653589793238462643383279502884;
    static constexpr double TWO_PI = 2.0 * PI;
    static constexpr double DEG_TO_RAD = PI / 180.0;
    static constexpr double RAD_TO_DEG = 180.0 / PI;
    
    /// Parametro gravitazionale Sole [AU³/day²] (OrbFit high precision)
    static constexpr double GM_SUN = 0.0002959122082855911;
    
    /// Obliquità eclittica J2000 [rad] (23° 26' 21".406)
    static constexpr double OBLIQUITY_J2000 = 0.4090926227129525;
    
    // ========================================================================
    // Conversioni principali
    // ========================================================================
    
    /**
     * @brief Converte elementi equinoziali in kepleriani
     * @param eq Elementi equinoziali (formato AstDyS)
     * @return Elementi kepleriani
     * 
     * Formule di conversione:
     * - e = sqrt(h² + k²)
     * - i = 2·atan(sqrt(p² + q²))
     * - Ω = atan2(p, q)
     * - ϖ = atan2(h, k)  [longitudine perielio]
     * - ω = ϖ - Ω
     * - M = λ - ϖ
     */
    static KeplerianElements equinoctialToKeplerian(
        const EquinoctialElements& eq);

    /**
     * @brief Converte elementi medi (AstDyS) in osculatori
     * @param mean Elementi medi kepleriani
     * @return Elementi osculatori kepleriani
     * 
     * Include:
     * - Perturbazioni periodiche planetarie (Milani-Knezevic)
     * - Correzioni a breve periodo (J2)
     */
    static KeplerianElements meanToOsculating(
        const KeplerianElements& mean);
    
    /**
     * @brief Converte elementi kepleriani in stato cartesiano (eclittico)
     * @param kep Elementi kepleriani
     * @return Stato cartesiano in frame eclittico J2000
     * 
     * Procedura:
     * 1. Risolve equazione Keplero per anomalia eccentrica E
     * 2. Calcola anomalia vera ν
     * 3. Calcola posizione e velocità nel piano orbitale
     * 4. Applica matrice di Gauss per rotazione in eclittico
     */
    static CartesianState keplerianToCartesian(
        const KeplerianElements& kep);
    
    /**
     * @brief Converte stato eclittico in ICRF (equatoriale J2000)
     * @param ecliptic Stato in coordinate eclittiche J2000
     * @return Stato in frame ICRF (equatoriale J2000)
     * 
     * Applica rotazione con obliquità ε = 23.439291°
     * Frame ICRF è compatibile con JPL Horizons
     */
    static CartesianState eclipticToICRF(
        const CartesianState& ecliptic);
    
    /**
     * @brief Converte ICRF in eclittico (operazione inversa)
     * @param icrf Stato in frame ICRF
     * @return Stato in coordinate eclittiche J2000
     */
    static CartesianState icrfToEcliptic(
        const CartesianState& icrf);
    
    // ========================================================================
    // Utility matematiche
    // ========================================================================
    
    /**
     * @brief Risolve equazione di Keplero: M = E - e·sin(E)
     * @param M Anomalia media [rad]
     * @param e Eccentricità
     * @param tol Tolleranza convergenza (default 1e-14)
     * @param max_iter Massimo iterazioni (default 20)
     * @return Anomalia eccentrica E [rad]
     * @throws std::runtime_error Se non converge
     * 
     * Metodo di Newton-Raphson:
     * E_{n+1} = E_n - (E_n - e·sin(E_n) - M) / (1 - e·cos(E_n))
     */
    static double solveKeplerEquation(double M, double e, 
                                     double tol = 1e-14,
                                     int max_iter = 20);
    
    /**
     * @brief Normalizza angolo in [0, 2π)
     * @param angle Angolo [rad]
     * @return Angolo normalizzato [rad]
     */
    static double normalizeAngle(double angle);
    
    /**
     * @brief Normalizza angolo in [-π, π)
     * @param angle Angolo [rad]
     * @return Angolo normalizzato [rad]
     */
    static double normalizeAngleSigned(double angle);
    
    // ========================================================================
    // Validazione
    // ========================================================================
    
    /**
     * @brief Valida stato cartesiano in frame ICRF
     * @param state Stato da validare
     * @return true se valido
     * 
     * Criteri:
     * - Coordinate finite
     * - Distanza in range [0.1, 100] AU
     * - Velocità < 100 AU/day
     */
    static bool validateICRF(const CartesianState& state);
    
    /**
     * @brief Calcola elementi kepleriani da stato cartesiano
     * @param state Stato cartesiano (ICRF o eclittico)
     * @param mu Parametro gravitazionale (default: GM_Sun)
     * @return Elementi kepleriani
     * 
     * Conversione inversa (stato → elementi)
     */
    static KeplerianElements cartesianToKeplerian(
        const CartesianState& state,
        double mu = GM_SUN);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_ORBITAL_CONVERSIONS_H
