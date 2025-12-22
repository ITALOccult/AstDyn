/**
 * @file phase1_candidate_screening.h
 * @brief Classe per FASE 1: Screening stelle candidate per occultazioni
 * 
 * Questa classe implementa la prima fase della ricerca di occultazioni:
 * 1. Crea un path ad alta risoluzione dell'asteroide
 * 2. Cerca stelle nel corridor usando UnifiedGaiaCatalog
 * 3. Calcola closest approach per ogni stella
 * 4. Filtra e restituisce solo le candidate interessanti
 * 
 * La Fase 2 (calcolo geometria occultazione precisa) usa queste candidate
 * come input per propagazioni ultra-precise.
 * 
 * @author IOccultCalc Team
 * @date 4 Dicembre 2025
 */

#ifndef PHASE1_CANDIDATE_SCREENING_H
#define PHASE1_CANDIDATE_SCREENING_H

#include <vector>
#include <memory>
#include <string>
#include <Eigen/Dense>

// Forward declarations FUORI dal namespace ioccultcalc
namespace astdyn {
    namespace propagation {
        class Propagator;
        class KeplerianElements;
    }
    namespace ephemeris {
        class PlanetaryEphemeris;
    }
}

namespace ioc {
    namespace gaia {
        class UnifiedGaiaCatalog;
    }
}

namespace ioccultcalc {

/**
 * @struct PathPoint
 * @brief Punto del path dell'asteroide
 */
struct PathPoint {
    double mjd_tdb;                 ///< Modified Julian Date (TDB)
    double ra_deg;                  ///< Right Ascension [degrees]
    double dec_deg;                 ///< Declination [degrees]
    Eigen::Vector3d pos_geo_au;     ///< Posizione geocentrica [AU, ICRF J2000]
    double distance_earth_au;       ///< Distanza dalla Terra [AU]
};

/**
 * @struct CandidateStar
 * @brief Stella candidata per occultazione (output Fase 1)
 */
struct CandidateStar {
    uint64_t source_id;                    ///< Gaia DR3 Source ID
    double ra_deg;                         ///< Right Ascension [degrees]
    double dec_deg;                        ///< Declination [degrees]
    double pm_ra_mas_yr;                   ///< Proper Motion RA [mas/yr]
    double pm_dec_mas_yr;                  ///< Proper Motion Dec [mas/yr]
    double phot_g_mean_mag;                ///< Gaia G magnitude
    double closest_approach_arcsec;        ///< Distanza minima dal path [arcsec]
    double closest_approach_mjd;           ///< Epoca del closest approach [MJD TDB]
    int closest_segment_index;             ///< Indice segmento path con CA minimo
    double angular_velocity_arcsec_per_sec; ///< Velocità angolare al CA [arcsec/sec]
};

/**
 * @struct Phase1Config
 * @brief Configurazione per lo screening Fase 1
 */
struct Phase1Config {
    // === PARAMETRI PATH ===
    double start_mjd_tdb;              ///< Inizio intervallo [MJD TDB]
    double end_mjd_tdb;                ///< Fine intervallo [MJD TDB]
    int path_interval_seconds;         ///< Intervallo tra punti path [secondi]
    
    // === PARAMETRI CORRIDOR ===
    double corridor_width_deg;         ///< Semi-larghezza corridor [gradi]
    double max_magnitude;              ///< Magnitudine limite
    double min_parallax;               ///< Parallasse minima [mas], -1 = no limit
    
    // === PARAMETRI FILTRO ===
    double closest_approach_threshold_arcsec; ///< Soglia CA per filtrare candidate [arcsec]
    
    // === PARAMETRI CATALOGO ===
    std::string catalog_mode;          ///< "online" o "local"
    
    // Costruttore con valori di default ottimizzati
    Phase1Config() 
        : start_mjd_tdb(0.0)
        , end_mjd_tdb(0.0)
        , path_interval_seconds(30)           // 30 sec = buon compromesso velocità/precisione
        , corridor_width_deg(0.0083)          // 30 arcsec = ~0.5 arcmin (RACCOMANDATO)
        , max_magnitude(18.0)                 // Stelle visibili da terra con setup amatoriale
        , min_parallax(-1.0)                  // Nessun limite su parallasse
        , closest_approach_threshold_arcsec(15.0)  // Solo CA < 15 arcsec (3x ombra tipica)
        , catalog_mode("local")               // Default: Catalogo locale multifile
    {}
    
    /**
     * @brief Configurazione conservativa (più stelle, più sicuro, più lento)
     */
    static Phase1Config conservative() {
        Phase1Config cfg;
        cfg.path_interval_seconds = 20;         // Path più denso
        cfg.corridor_width_deg = 0.0167;        // 60 arcsec = 1 arcmin
        cfg.closest_approach_threshold_arcsec = 30.0;
        return cfg;
    }
    
    /**
     * @brief Configurazione veloce (meno stelle, più veloce, rischio di perdere candidate)
     */
    static Phase1Config fast() {
        Phase1Config cfg;
        cfg.path_interval_seconds = 60;         // Path meno denso
        cfg.corridor_width_deg = 0.0056;        // 20 arcsec
        cfg.closest_approach_threshold_arcsec = 10.0;
        return cfg;
    }
    
    /**
     * @brief Configurazione per survey completo (trova tutto, molto lento)
     */
    static Phase1Config survey() {
        Phase1Config cfg;
        cfg.path_interval_seconds = 10;         // Path molto denso
        cfg.corridor_width_deg = 0.05;          // 180 arcsec = 3 arcmin
        cfg.closest_approach_threshold_arcsec = 60.0;
        return cfg;
    }
};

/**
 * @struct Phase1Results
 * @brief Risultati della Fase 1
 */
struct Phase1Results {
    std::vector<PathPoint> path;                ///< Path completo dell'asteroide
    std::vector<CandidateStar> all_stars;       ///< Tutte le stelle nel corridor
    std::vector<CandidateStar> candidates;      ///< Stelle candidate filtrate (CA < threshold)
    
    // Statistiche
    double propagation_time_ms;                 ///< Tempo propagazione path [ms]
    double corridor_query_time_ms;              ///< Tempo query catalogo [ms]
    double closest_approach_calc_time_ms;       ///< Tempo calcolo CA [ms]
    int num_path_points;                        ///< Numero punti nel path
    int num_stars_in_corridor;                  ///< Numero stelle nel corridor
    int num_candidates_filtered;                ///< Numero candidate filtrate
    
    Phase1Results()
        : propagation_time_ms(0), corridor_query_time_ms(0)
        , closest_approach_calc_time_ms(0), num_path_points(0)
        , num_stars_in_corridor(0), num_candidates_filtered(0)
    {}
};

/**
 * @class Phase1CandidateScreening
 * @brief Implementa lo screening di Fase 1 per ricerca occultazioni
 * 
 * Workflow tipico:
 * @code
 * // 1. Configura
 * Phase1Config config = Phase1Config::conservative();
 * config.start_mjd_tdb = 61007.0;
 * config.end_mjd_tdb = 61008.0;
 * 
 * // 2. Crea screener
 * Phase1CandidateScreening screener;
 * 
 * // 3. Carica elementi orbitali
 * screener.loadAsteroidFromEQ1("17030_astdys.eq1");
 * 
 * // 4. Esegui screening
 * auto results = screener.screenCandidates(config);
 * 
 * // 5. Usa candidate per Fase 2
 * for (const auto& candidate : results.candidates) {
 *     // Fase 2: calcolo geometria occultazione precisa
 *     computeOccultationGeometry(candidate);
 * }
 * @endcode
 */
class Phase1CandidateScreening {
public:
    Phase1CandidateScreening();
    ~Phase1CandidateScreening();
    
    // Non copiabile
    Phase1CandidateScreening(const Phase1CandidateScreening&) = delete;
    Phase1CandidateScreening& operator=(const Phase1CandidateScreening&) = delete;
    
    /**
     * @brief Carica elementi orbitali dell'asteroide da file .eq1
     * @param eq1_path Path al file OrbFit .eq1
     * @return true se caricamento riuscito
     */
    bool loadAsteroidFromEQ1(const std::string& eq1_path);
    
    /**
     * @brief Carica elementi orbitali dal database JSON locale
     * @param asteroid_number Numero dell'asteroide (es: 17030)
     * @param json_path Path al file JSON (default: ~/.ioccultcalc/data/all_numbered_asteroids.json)
     * @return true se caricamento riuscito
     */
    bool loadAsteroidFromJSON(int asteroid_number, const std::string& json_path = "");
    
    /**
     * @brief Imposta elementi orbitali direttamente
     * @param elements Elementi kepleriani
     */
    void setOrbitalElements(const astdyn::propagation::KeplerianElements& elements);
    
    /**
     * @brief Esegue lo screening completo di Fase 1
     * @param config Configurazione dello screening
     * @return Risultati con path, stelle, e candidate filtrate
     * @throw std::runtime_error se elementi non caricati o catalogo non inizializzato
     */
    Phase1Results screenCandidates(const Phase1Config& config);
    
    /**
     * @brief Ottiene gli elementi orbitali correnti
     * @return Elementi kepleriani caricati
     */
    const astdyn::propagation::KeplerianElements& getOrbitalElements() const;
    
    /**
     * @brief Verifica se elementi orbitali sono stati caricati
     */
    bool hasOrbitalElements() const;
    
    /**
     * @brief Imposta il catalogo Gaia da usare (opzionale, usa getInstance se non impostato)
     * @param catalog Puntatore al catalogo
     */
    void setCatalog(ioc::gaia::UnifiedGaiaCatalog* catalog);

    /**
     * @brief Imposta il livello di verbosità
     * @param level 0=quiet, 1=minimal, 2=full
     */
    void setVerbose(int level);

    
private:
    class Impl;
    std::unique_ptr<Impl> pimpl_;
    
    /**
     * @brief Crea path ad alta risoluzione
     */
    std::vector<PathPoint> createHighResolutionPath(const Phase1Config& config, double& time_ms);
    
    /**
     * @brief Query corridor nel catalogo Gaia
     */
    std::vector<CandidateStar> queryCorridor(const std::vector<PathPoint>& path,
                                              const Phase1Config& config,
                                              double& time_ms);
    
    /**
     * @brief Calcola closest approach per ogni stella
     */
    void computeClosestApproaches(const std::vector<PathPoint>& path,
                                   std::vector<CandidateStar>& stars,
                                   double& time_ms);
    
    /**
     * @brief Filtra stelle per closest approach
     */
    std::vector<CandidateStar> filterCandidates(const std::vector<CandidateStar>& stars,
                                                 double threshold_arcsec);
    
    /**
     * @brief Calcola velocità angolare al closest approach
     */
    double computeAngularVelocity(const std::vector<PathPoint>& path,
                                   const CandidateStar& star);
};

} // namespace ioccultcalc

#endif // PHASE1_CANDIDATE_SCREENING_H
