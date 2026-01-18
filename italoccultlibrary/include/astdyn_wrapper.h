/**
 * @file astdyn_wrapper.h
 * @brief Wrapper semplificato per AstDyn library
 * @author IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 * 
 * Wrapper che incapsula l'uso di AstDyn per la propagazione orbitale.
 * Basato sull'API reale di AstDyn v1.0.0
 */

#ifndef IOCCULTCALC_ASTDYN_WRAPPER_H
#define IOCCULTCALC_ASTDYN_WRAPPER_H

#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/io/parsers/OrbFitEQ1Parser.hpp>
#include <astdyn/api/OrbitFitAPI.hpp>
#include <astdyn/propagation/HighPrecisionPropagator.hpp>
#include <astdyn/propagation/STMPropagator.hpp>
#include "ioccultcalc/orbital_elements.h"
#include <astdyn/core/Types.hpp>
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>
#include <optional>

namespace ioccultcalc {

/**
 * @struct PropagationSettings
 * @brief Configurazione per la propagazione orbitale
 */
struct PropagationSettings {
    double tolerance = 1e-12;        ///< Tolleranza RKF78 [AU]
    double initial_step = 0.1;       ///< Step iniziale [giorni]
    bool include_planets = true;     ///< Perturbazioni planetarie
    bool include_relativity = true;  ///< Relatività generale
    bool include_asteroids = true;   ///< Perturbazioni asteroidi
    bool perturb_mercury = true;
    bool perturb_venus = true;
    bool perturb_earth = true;
    bool perturb_mars = true;
    bool perturb_jupiter = true;
    bool perturb_saturn = true;
    bool perturb_uranus = true;
    bool perturb_neptune = true;
    
    /// Preset: massima accuratezza (JPL-compliant)
    static PropagationSettings highAccuracy() {
        PropagationSettings s;
        s.tolerance = 1e-12;
        s.initial_step = 0.1;
        return s;  // Tutti i flag già true di default
    }
    
    /// Preset: veloce (screening)
    static PropagationSettings fast() {
        PropagationSettings s;
        s.tolerance = 1e-9;
        s.initial_step = 0.5;
        s.include_relativity = false;
        s.include_asteroids = false;
        s.perturb_uranus = false;
        s.perturb_neptune = false;
        return s;
    }
};

/**
 * @struct CartesianStateICRF
 * @brief Stato cartesiano in frame ICRF (equatoriale J2000)
 */
struct CartesianStateICRF {
    Eigen::Vector3d position;  ///< Posizione [AU]
    Eigen::Vector3d velocity;  ///< Velocità [AU/day]
    double epoch_mjd_tdb;      ///< Epoca [MJD TDB]
    std::optional<Eigen::Matrix<double, 6, 6>> covariance; ///< Covarianza [6x6]
    
    CartesianStateICRF() : epoch_mjd_tdb(0.0) {
        position = Eigen::Vector3d::Zero();
        velocity = Eigen::Vector3d::Zero();
    }
};

/**
 * @struct FitResult
 * @brief Risultato del fitting orbitale (Orbit Determination)
 */
struct FitResult {
    bool success = false;
    std::string message;
    double rms_ra = 0.0;     ///< [arcsec]
    double rms_dec = 0.0;    ///< [arcsec]
    int iterations = 0;
    int num_observations = 0;
    bool converged = false;
};

/**
 * @class AstDynWrapper
 * @brief Wrapper semplificato per propagazione con AstDyn
 * 
 * Esempio d'uso:
 * @code
 *   // Crea wrapper con configurazione default
 *   AstDynWrapper wrapper;
 *   
 *   // Carica elementi da file .eq1
 *   wrapper.loadFromEQ1File("17030.eq1");
 *   
 *   // Propaga a epoch target
 *   auto state = wrapper.propagateToEpoch(61007.0);  // MJD TDB
 *   
 *   std::cout << "Posizione: " << state.position.transpose() << std::endl;
 * @endcode
 */
class AstDynWrapper {
public:
    /**
     * @brief Costruttore con configurazione
     * @param settings Impostazioni propagazione (default: highAccuracy)
     */
    explicit AstDynWrapper(const PropagationSettings& settings = PropagationSettings::highAccuracy());
    
    /**
     * @brief Distruttore
     */
    ~AstDynWrapper();
    
    /**
     * @brief Carica elementi orbitali da file .eq1
     * @param filepath Percorso file .eq1 (formato OEF2.0 AstDyS)
     * @return true se caricamento riuscito
     */
    bool loadFromEQ1File(const std::string& filepath);
    
    /**
     * @brief Inietta manualmente elementi orbitali (JPL Injector)
     * 
     * @param a Semiasse maggiore [AU]
     * @param e Eccentricità
     * @param i Inclinazione [rad]
     * @param Omega Longitudine nodo ascendente [rad]
     * @param omega Argomento perielio [rad]
     * @param M Anomalia media [rad]
     * @param epoch Epoca [MJD TDB]
     * @param name Nome oggetto
     */
    void injectElements(double a, double e, double i, double Omega, double omega, double M, double epoch, const std::string& name = "INJECTED", double H=0, double G=0.15, double dia=0) {
        setKeplerianElements(a, e, i, Omega, omega, M, epoch, name, 
                             astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC,
                             H, G, dia);
    }

    /**
     * @brief Carica elementi orbitali osculanti da JPL Horizons
     * 
     * @param designation Nome o numero asteroide (es. "13477")
     * @param epoch_mjd_tdb Epoca richiesta [MJD TDB]
     * @return true se operazione riuscita (richiede connessione internet)
     */
    bool loadFromHorizons(const std::string& designation, double epoch_mjd_tdb);

    /**
     * @brief Imposta elementi kepleriani manualmente
     * @param a Semiasse maggiore [AU]
     * @param e Eccentricità
     * @param i Inclinazione [rad]
     * @param Omega Longitudine nodo ascendente [rad]
     * @param omega Argomento perielio [rad]
     * @param M Anomalia media [rad]
     * @param epoch_mjd_tdb Epoca [MJD TDB]
     * @param name Nome oggetto (opzionale)
     * @param frame Frame di riferimento degli elementi (default: ECLIPTIC)
     */
    void setKeplerianElements(double a, double e, double i, 
                              double Omega, double omega, double M,
                              double epoch_mjd_tdb,
                              const std::string& name = "",
                              astdyn::propagation::HighPrecisionPropagator::InputFrame frame = 
                                  astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC,
                              double H = 0.0, double G = 0.15, double dia = 0.0,
                              const std::optional<Eigen::Matrix<double, 6, 6>>& covariance = std::nullopt);
    
    /**
     * @brief Imposta elementi cartesiani manualmente
     * @param pos Posizione [AU]
     * @param vel Velocità [AU/day]
     * @param epoch_mjd_tdb Epoca [MJD TDB]
     * @param name Nome oggetto
     * @param frame Frame di riferimento (default: EQUATORIAL/ICRF)
     */
    void setCartesianElements(const Eigen::Vector3d& pos, const Eigen::Vector3d& vel,
                             double epoch_mjd_tdb, const std::string& name = "",
                             astdyn::propagation::HighPrecisionPropagator::InputFrame frame = 
                                 astdyn::propagation::HighPrecisionPropagator::InputFrame::EQUATORIAL);

    /**
     * @brief Imposta gli elementi orbitali da struct ioccultcalc.
     * Applica le conversioni necessarie (Mean -> Osculating) se il tipo è MEAN_ASTDYS.
     */
    void setAsteroidElements(const AstDynEquinoctialElements& elements);
    
    /**
     * @brief Imposta elementi orbitali equinoziali (OSCULANTI)
     * @param a Semiasse maggiore [AU]
     * @param h e*sin(omega+Omega)
     * @param k e*cos(omega+Omega)
     * @param p tan(i/2)*sin(Omega)
     * @param q tan(i/2)*cos(Omega)
     * @param lambda Longitudine media [rad]
     * @param epoch_mjd_tdb Epoca [MJD TDB]
     * @param name Nome oggetto (opzionale)
     * @param frame Frame di riferimento (default: ECLIPTIC, assumendo input osculanti)
     */
    void setEquinoctialElements(double a, double h, double k,
                                double p, double q, double lambda,
                                double epoch_mjd_tdb,
                                const std::string& name = "",
                                astdyn::propagation::HighPrecisionPropagator::InputFrame frame = 
                                    astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC,
                                const std::optional<Eigen::Matrix<double, 6, 6>>& covariance = std::nullopt);

    /**
     * @brief Imposta elementi equinoziali MEDI (es. da allnum.cat / AstDyS)
     * 
     * Esegue automaticamente la conversione Mean -> Osculating (ICRF)
     * richiesta per l'alta precisione.
     * 
     * @param a Semiasse maggiore [AU]
     * @param h e*sin(omega+Omega)
     * @param k e*cos(omega+Omega)
     * @param p tan(i/2)*sin(Omega)
     * @param q tan(i/2)*cos(Omega)
     * @param lambda Longitudine media [rad]
     * @param epoch_mjd_tdb Epoca [MJD TDB]
     * @param name Nome oggetto (opzionale)
     */
    void setMeanEquinoctialElements(double a, double h, double k,
                                    double p, double q, double lambda,
                                    double epoch_mjd_tdb,
                                    const std::string& name = "");
    
    /**
     * @brief Propaga orbita a epoca target
     * @param target_mjd_tdb Epoca target [MJD TDB]
     * @return Stato cartesiano in frame ICRF
     * @throws std::runtime_error Se elementi non inizializzati
     */
    CartesianStateICRF propagateToEpoch(double target_mjd_tdb);

    /**
     * @brief Esegue un fitting orbitale (Least Squares)
     * @param rwo_filepath Percorso file osservazioni (.rwo)
     * @param verbose Abilita output dettagliato
     * @return Risultato del fit
     */
    FitResult runFit(const std::string& rwo_filepath, bool verbose = true);

    /**
     * @brief Calcola osservazione geocentrica (RA/Dec) ad alta precisione
     * @param target_mjd_tdb Epoca target [MJD TDB]
     * @return Risultato osservazione (gradi, AU, ecc.)
     */
    astdyn::propagation::HighPrecisionPropagator::ObservationResult 
    calculateObservation(double target_mjd_tdb);
    
    /**
     * @brief Ottiene epoca corrente elementi orbitali
     * @return Epoca [MJD TDB]
     */
    double getCurrentEpoch() const { return current_epoch_mjd_; }
    
    std::string getObjectName() const;
    
    /**
     * @brief Verifica se elementi sono stati caricati
     * @return true se elementi validi disponibili
     */
    bool isInitialized() const { return initialized_; }
    
    /**
     * @brief Ottiene statistiche ultima propagazione
     * @return Stringa con informazioni (step, tempo, ecc.)
     */
    std::string getLastPropagationStats() const { return last_stats_; }

    /**
     * @struct ApparentState
     * @brief Stato apparente geocentrico per compatibilità Phase2
     */
    struct ApparentState {
        double ra_deg;
        double dec_deg;
        double distance_au;
        Eigen::Vector3d position;
    };

    /**
     * @brief Ottiene stato apparente geocentrico (RA/Dec) a epoca MJD
     * @param mjd_tdb Epoca [MJD TDB]
     * @return Stato apparente
     */
    ApparentState getApparentStateGeocentric(double mjd_tdb);

    /**
     * @brief Ottiene elementi kepleriani correnti (incluso Mag)
     * @return Elementi kepleriani (IOrbitParser format)
     */
    astdyn::io::IOrbitParser::OrbitalElements getKeplerianElements() const;

    // Compatibilità legacy
    void setSPKReader(std::shared_ptr<void> reader);
    
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> getPlanetaryEphemeris() const { return ephemeris_; }

private:
    PropagationSettings settings_;
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephemeris_;
    std::unique_ptr<astdyn::propagation::Propagator> propagator_;
    std::unique_ptr<astdyn::propagation::HighPrecisionPropagator> high_prop_;
    
    astdyn::io::IOrbitParser::OrbitalElements current_elements_;
    astdyn::propagation::HighPrecisionPropagator::InputFrame current_frame_;
    std::optional<Eigen::Matrix<double, 6, 6>> current_covariance_;
    double current_epoch_mjd_;
    std::string object_name_;
    bool initialized_;
    std::string last_stats_;
    
    /// Inizializza propagatore con settings correnti
    void initializePropagator();
    
    static std::string extractObjectNameFromEQ1(const std::string& filepath);

    /// Conversione frame ECLM J2000 → ICRF
    static CartesianStateICRF eclipticToICRF(const astdyn::propagation::CartesianElements& ecl_state,
                                             double epoch_mjd_tdb);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_ASTDYN_WRAPPER_H
