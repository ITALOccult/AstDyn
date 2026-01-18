/**
 * @file astdyn_wrapper.cpp
 * @brief Implementazione del wrapper AstDyn
 * @author IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 */

#include "astdyn_wrapper.h"
#include "orbital_conversions.h"
#include "ioccultcalc/jpl_horizons_client.h"
#include <sstream>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

namespace ioccultcalc {

// Helper per estrarre il nome dell'oggetto dal file .eq1
static std::string extractObjectNameFromEQ1(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        return "";
    }
    
    bool header_ended = false;
    std::string line;
    while (std::getline(file, line)) {
        // Cerca fine dell'header
        if (line.find("END_OF_HEADER") != std::string::npos) {
            header_ended = true;
            continue;
        }
        
        // Dopo END_OF_HEADER, la prima linea non-commento è il nome/numero
        if (header_ended) {
            // Salta righe vuote e commenti
            if (line.empty() || line[0] == '!') {
                continue;
            }
            
            // Rimuovi spazi bianchi iniziali e finali
            size_t start = line.find_first_not_of(" \t\r\n");
            size_t end = line.find_last_not_of(" \t\r\n");
            
            if (start != std::string::npos && end != std::string::npos) {
                return line.substr(start, end - start + 1);
            }
        }
    }
    
    return "";
}

AstDynWrapper::AstDynWrapper(const PropagationSettings& settings)
    : settings_(settings)
    , current_epoch_mjd_(0.0)
    , initialized_(false)
    , current_frame_(astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC)
{
    std::cout << "[AstDynWrapper] Constructor started." << std::endl;
    // Crea effemeridi planetarie (usa DE440 embedded di default)
    std::cout << "  - Creating PlanetaryEphemeris..." << std::endl;
    ephemeris_ = std::make_shared<astdyn::ephemeris::PlanetaryEphemeris>();
    
    // Inizializza propagatore
    std::cout << "  - Initializing propagator..." << std::endl;
    initializePropagator();
    std::cout << "[AstDynWrapper] Constructor finished." << std::endl;
}

AstDynWrapper::~AstDynWrapper() = default;

void AstDynWrapper::initializePropagator() {
    std::cout << "  - Creating RKF78Integrator..." << std::endl;
    // Crea integratore RKF78
    auto integrator = std::make_unique<astdyn::propagation::RKF78Integrator>(
        settings_.initial_step,
        settings_.tolerance
    );
    
    // Configura perturbazioni
    astdyn::propagation::PropagatorSettings prop_settings;
    prop_settings.include_planets = settings_.include_planets;
    prop_settings.include_relativity = settings_.include_relativity;
    prop_settings.include_asteroids = settings_.include_asteroids;
    prop_settings.perturb_mercury = settings_.perturb_mercury;
    prop_settings.perturb_venus = settings_.perturb_venus;
    prop_settings.perturb_earth = settings_.perturb_earth;
    prop_settings.perturb_mars = settings_.perturb_mars;
    prop_settings.perturb_jupiter = settings_.perturb_jupiter;
    prop_settings.perturb_saturn = settings_.perturb_saturn;
    prop_settings.perturb_uranus = settings_.perturb_uranus;
    prop_settings.perturb_neptune = settings_.perturb_neptune;
    
    // Crea propagatore
    std::cout << "  - Creating Propagator object..." << std::endl;
    propagator_ = std::make_unique<astdyn::propagation::Propagator>(
        std::move(integrator),
        ephemeris_,
        prop_settings
    );
    std::cout << "  - Propagator object created." << std::endl;
}

bool AstDynWrapper::loadFromEQ1File(const std::string& filepath) {
    try {
        // 1. Usa OrbitFitAPI per il solo parsing degli elementi equinoziali (Medi, Eclittici)
        auto equ = astdyn::api::OrbitFitAPI::parse_eq1(filepath);
        
        // 2. Converti Equinoziali -> Kepleriani (sempre Medi, Eclittici)
        auto kep_mean = astdyn::propagation::equinoctial_to_keplerian(equ);
        
        // 3. Converti Medi -> Osculanti (Milani-Knezevic)
        // Nota: mean_to_osculating applica le correzioni MK se a è nel range main belt
        auto kep_osc = astdyn::propagation::mean_to_osculating(kep_mean);
        
        std::cout << "[AstDynWrapper] Computed Osculating Elements (Epoch " << kep_osc.epoch_mjd_tdb << "):\n"
                  << "  a=" << kep_osc.semi_major_axis << " e=" << kep_osc.eccentricity << " i=" << kep_osc.inclination << "\n"
                  << "  Omega=" << kep_osc.longitude_ascending_node << " omega=" << kep_osc.argument_perihelion << " M=" << kep_osc.mean_anomaly << "\n";
        
        // Salva elementi internamente (Osculanti, Eclittici)
        current_elements_.semi_major_axis = kep_osc.semi_major_axis;
        current_elements_.eccentricity = kep_osc.eccentricity;
        current_elements_.inclination = kep_osc.inclination;
        current_elements_.longitude_asc_node = kep_osc.longitude_ascending_node;
        current_elements_.argument_perihelion = kep_osc.argument_perihelion;
        current_elements_.mean_anomaly = kep_osc.mean_anomaly;
        current_elements_.epoch_mjd_tdb = kep_osc.epoch_mjd_tdb;
        current_elements_.gravitational_parameter = kep_osc.gravitational_parameter;
        
        current_epoch_mjd_ = kep_osc.epoch_mjd_tdb;
        
        // IMPORTANTE: Utilizziamo ECLIPTIC per coerenza con Phase 1 
        // e per evitare bug di rotazione in OrbitFitAPI.
        // HighPrecisionPropagator ruoterà internamente in Equatorial (ICRF).
        current_frame_ = astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC;
        
        // Il parser non estrae object_name, lo leggiamo manualmente
        object_name_ = extractObjectNameFromEQ1(filepath);
        if (object_name_.empty()) object_name_ = "Unknown";
        
        std::cout << "[AstDynWrapper:" << this << "] loadFromEQ1File: asteroid=" << object_name_ 
                  << " (Frame: ECLIPTIC, Mean->Osc: YES) setting initialized_ = true" << std::endl;
        initialized_ = true;
        // DEFERRED: initializePropagator(); // Avoid slow initialization here
        
        // Nota: Poiché abbiamo convertito in Equatorial (ICRF), 
        // dobbiamo assicurarci che il propagatore sia configurato per EQUATORIAL
        // in propagateToEpoch (implementato via propagate_keplerian).
        
        return true;
        
    } catch (const std::exception& e) {
        std::cout << "[AstDynWrapper:" << this << "] loadFromEQ1File: FAILED, setting initialized_ = false" << std::endl;
        initialized_ = false;
        last_stats_ = std::string("Errore caricamento alta precisione (convert_mean_to_osculating): ") + e.what();
        return false;
    }
}

bool AstDynWrapper::loadFromHorizons(const std::string& designation, double epoch_mjd_tdb) {
    try {
        JPLHorizonsClient client;
        JulianDate epoch;
        epoch.jd = epoch_mjd_tdb + 2400000.5;
        
        std::cout << "[AstDynWrapper] Fetching elements from JPL Horizons for " << designation << " at MJD " << epoch_mjd_tdb << std::endl;
        
        auto elem = client.getOsculatingElements(designation, epoch);
        
        // Imposta gli elementi (Horizons restituisce elementi OSCULANTI in frame ECLIPTIC J2000 per default se non specificato altrimenti)
        // La nostra JPLHorizonsClient::getOsculatingElements imposta REF_PLANE='ECLIPTIC'
        setKeplerianElements(elem.a, elem.e, elem.i, elem.Omega, elem.omega, elem.M, 
                             elem.epoch.jd - 2400000.5, designation, 
                             astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC);
        
        std::cout << "[AstDynWrapper] loadFromHorizons: SUCCESS" << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "[AstDynWrapper] loadFromHorizons failed: " << e.what() << std::endl;
        return false;
    }
}

void AstDynWrapper::setKeplerianElements(
    double a, double e, double i,
    double Omega, double omega, double M,
    double epoch_mjd_tdb,
    const std::string& name,
    astdyn::propagation::HighPrecisionPropagator::InputFrame frame,
    double H, double G, double dia,
    const std::optional<Eigen::Matrix<double, 6, 6>>& covariance)
{
    // Salva elementi per uso futuro
    current_elements_.semi_major_axis = a;
    current_elements_.eccentricity = e;
    current_elements_.inclination = i;
    current_elements_.longitude_asc_node = Omega;
    current_elements_.argument_perihelion = omega;
    current_elements_.mean_anomaly = M;
    current_elements_.epoch_mjd_tdb = epoch_mjd_tdb;
    current_elements_.object_name = name;
    current_elements_.magnitude = H;
    current_elements_.mag_slope = G;
    current_elements_.diameter = dia;
    current_elements_.gravitational_parameter = astdyn::constants::GMS; // Default if not provided
    
    current_epoch_mjd_ = epoch_mjd_tdb;
    current_frame_ = frame;
    object_name_ = name;
    current_covariance_ = covariance;
    std::cout << "[AstDynWrapper:" << this << "] setKeplerianElements: asteroid=" << name << " state set to INITIALIZED (Covariance: " << (covariance.has_value()?"YES":"NO") << ")." << std::endl;
    std::cout << "[AstDynWrapper]   Processing Frame: " << (frame == astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC ? "ECLIPTIC" : "EQUATORIAL") << "\n";
    initialized_ = true;
    
    // Re-inizializza il propagatore se necessario (opzionale, propagateToEpoch lo usa)
    initializePropagator();
}

void AstDynWrapper::setAsteroidElements(const AstDynEquinoctialElements& elements) {
    if (elements.type == ElementType::MEAN_ASTDYS) {
        // Applica correzione MILANI-KNEZEVIC via OrbitalConversions
        auto osc = elements.toOsculatingKeplerian();
        
        // Assumiamo che se abbiamo covarianza sia già nel frame target (ECLIPTIC di solito per MEAN)
        std::optional<Eigen::Matrix<double, 6, 6>> cov_opt;
        if (elements.hasCovariance && elements.covariance.size() == 21) {
            Eigen::Matrix<double, 6, 6> cov = Eigen::Matrix<double, 6, 6>::Zero();
            int idx = 0;
            for (int i = 0; i < 6; ++i) {
                for (int j = i; j < 6; ++j) {
                    cov(i, j) = elements.covariance[idx++];
                    cov(j, i) = cov(i, j);
                }
            }
            cov_opt = cov;
        }

        setKeplerianElements(osc.a, osc.e, osc.i, osc.Omega, osc.omega, osc.M, 
                            osc.epoch.jd - 2400000.5, elements.name, 
                            (elements.frame == FrameType::EQUATORIAL_ICRF) ? 
                             astdyn::propagation::HighPrecisionPropagator::InputFrame::EQUATORIAL : 
                             astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC,
                             elements.H, elements.G, elements.diameter,
                             cov_opt);
        
        // Importante: toOsculatingKeplerian restituisce ElementType::OSCULATING
        // ma noi ricordiamo il frame originale se possibile.
    } else {
        // Già osculanti o ignoti, carica direttamente
        auto kep = elements.toKeplerian();
        auto target_frame = (elements.frame == FrameType::EQUATORIAL_ICRF) ? 
                              astdyn::propagation::HighPrecisionPropagator::InputFrame::EQUATORIAL : 
                              astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC;

        std::cout << "[AstDynWrapper] setAsteroidElements DEBUG: "
                  << "Name=" << elements.name 
                  << " FrameType=" << (int)elements.frame
                  << " -> InputFrame=" << (int)target_frame
                  << " omega=" << kep.omega << " Omega=" << kep.Omega << " M=" << kep.M << "\n";

        std::optional<Eigen::Matrix<double, 6, 6>> cov_opt;
        if (elements.hasCovariance && elements.covariance.size() == 21) {
            Eigen::Matrix<double, 6, 6> cov = Eigen::Matrix<double, 6, 6>::Zero();
            int idx = 0;
            for (int i = 0; i < 6; ++i) {
                for (int j = i; j < 6; ++j) {
                    cov(i, j) = elements.covariance[idx++];
                    cov(j, i) = cov(i, j);
                }
            }
            cov_opt = cov;
        }

        setKeplerianElements(kep.a, kep.e, kep.i, kep.Omega, kep.omega, kep.M, 
                            kep.epoch.jd - 2400000.5, elements.name,
                            target_frame,
                             elements.H, elements.G, elements.diameter,
                             cov_opt);
    }
}

void AstDynWrapper::setCartesianElements(const Eigen::Vector3d& pos, const Eigen::Vector3d& vel,
                                         double epoch_mjd_tdb, const std::string& name,
                                         astdyn::propagation::HighPrecisionPropagator::InputFrame frame)
{
    initialized_ = true;
    object_name_ = name;
    current_epoch_mjd_ = epoch_mjd_tdb;
    current_frame_ = frame;

    // Converti in elementi Kepleriani per consistenza con current_elements_
    astdyn::propagation::CartesianElements cart;
    cart.position = pos;
    cart.velocity = vel;
    cart.epoch_mjd_tdb = epoch_mjd_tdb;
    cart.gravitational_parameter = astdyn::constants::GMS;

    auto kep = astdyn::propagation::cartesian_to_keplerian(cart);
    
    current_elements_.object_name = name;
    current_elements_.epoch_mjd_tdb = epoch_mjd_tdb;
    current_elements_.semi_major_axis = kep.semi_major_axis;
    current_elements_.eccentricity = kep.eccentricity;
    current_elements_.inclination = kep.inclination;
    current_elements_.longitude_asc_node = kep.longitude_ascending_node;
    current_elements_.argument_perihelion = kep.argument_perihelion;
    current_elements_.mean_anomaly = kep.mean_anomaly;
    current_elements_.gravitational_parameter = astdyn::constants::GMS;

    std::cout << "[AstDynWrapper:" << this << "] setCartesianElements: asteroid=" << name 
              << " epoch=" << epoch_mjd_tdb << " frame=" << (int)frame << " state set to INITIALIZED.\n";
}

void AstDynWrapper::setEquinoctialElements(
    double a, double h, double k,
    double p, double q, double lambda,
    double epoch_mjd_tdb,
    const std::string& name,
    astdyn::propagation::HighPrecisionPropagator::InputFrame frame,
    const std::optional<Eigen::Matrix<double, 6, 6>>& covariance)
{
    // Converti equinoziali in kepleriani (mantiene il frame logico, conversione matematica pura)
    astdyn::propagation::EquinoctialElements eq;
    eq.a = a;
    eq.h = h;
    eq.k = k;
    eq.p = p;
    eq.q = q;
    eq.lambda = lambda;
    
    auto kep = astdyn::propagation::equinoctial_to_keplerian(eq);
    
    // Salva come kepleriani internamente
    current_elements_.semi_major_axis = kep.semi_major_axis;
    current_elements_.eccentricity = kep.eccentricity;
    current_elements_.inclination = kep.inclination;
    current_elements_.longitude_asc_node = kep.longitude_ascending_node;
    current_elements_.argument_perihelion = kep.argument_perihelion;
    current_elements_.mean_anomaly = kep.mean_anomaly;
    current_elements_.epoch_mjd_tdb = epoch_mjd_tdb;
    current_elements_.object_name = name;
    
    current_epoch_mjd_ = epoch_mjd_tdb;
    current_frame_ = frame;
    object_name_ = name;
    current_covariance_ = covariance;
    initialized_ = true;
    
    initializePropagator();
}

void AstDynWrapper::setMeanEquinoctialElements(
    double a, double h, double k,
    double p, double q, double lambda,
    double epoch_mjd_tdb,
    const std::string& name)
{
    // 1. Costruisci struttura EquinoctialElements per API (namespace propagation)
    astdyn::propagation::EquinoctialElements mean_elements;
    mean_elements.a = a;
    mean_elements.h = h;
    mean_elements.k = k;
    mean_elements.p = p;
    mean_elements.q = q;
    mean_elements.lambda = lambda;
    mean_elements.epoch_mjd_tdb = epoch_mjd_tdb;
    mean_elements.gravitational_parameter = astdyn::constants::GMS;
    
    // 2. Converti in Osculanti (ICRF) usando OrbitFitAPI
    auto kep_osc = astdyn::api::OrbitFitAPI::convert_mean_equinoctial_to_osculating(mean_elements);
    
    // 3. Salva elementi osculanti
    current_elements_.semi_major_axis = kep_osc.semi_major_axis;
    current_elements_.eccentricity = kep_osc.eccentricity;
    current_elements_.inclination = kep_osc.inclination;
    current_elements_.longitude_asc_node = kep_osc.longitude_ascending_node;
    current_elements_.argument_perihelion = kep_osc.argument_perihelion;
    current_elements_.mean_anomaly = kep_osc.mean_anomaly;
    current_elements_.epoch_mjd_tdb = kep_osc.epoch_mjd_tdb;
    current_elements_.object_name = name;
    
    current_epoch_mjd_ = kep_osc.epoch_mjd_tdb;
    
    // IMPORTANTE: Poiché convert_mean_equinoctial_to_osculating restituisce ICRF Equatoriali
    current_frame_ = astdyn::propagation::HighPrecisionPropagator::InputFrame::EQUATORIAL;
    
    object_name_ = name;
    std::cout << "[AstDynWrapper] setMeanEquinoctialElements: asteroid=" << name << " state set to INITIALIZED." << std::endl;
    initialized_ = true;
    
    initializePropagator();
}

CartesianStateICRF AstDynWrapper::propagateToEpoch(double target_mjd_tdb) {
    if (!initialized_) {
        throw std::runtime_error("AstDynWrapper: elementi non inizializzati. "
                                 "Chiamare loadFromEQ1File() o setKeplerianElements() prima.");
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::cout << "  - Preparing KeplerianElements..." << std::endl;
    // Crea KeplerianElements per propagazione (struct in propagation namespace)
    astdyn::propagation::KeplerianElements kep_initial;
    kep_initial.epoch_mjd_tdb = current_elements_.epoch_mjd_tdb;
    kep_initial.semi_major_axis = current_elements_.semi_major_axis;
    kep_initial.eccentricity = current_elements_.eccentricity;
    kep_initial.inclination = current_elements_.inclination;
    kep_initial.longitude_ascending_node = current_elements_.longitude_asc_node;
    kep_initial.argument_perihelion = current_elements_.argument_perihelion;
    kep_initial.mean_anomaly = current_elements_.mean_anomaly;
    kep_initial.gravitational_parameter = astdyn::constants::GMS;  // AU³/day²
    
    std::cout << "  - Propagating with AstDyn..." << std::endl;
    if (!propagator_) {
        std::cout << "  - ERROR: Propagator is NULL!" << std::endl;
        throw std::runtime_error("Propagator is NULL");
    }

    Eigen::Vector3d pos_final, vel_final;
    std::optional<Eigen::Matrix<double, 6, 6>> cov_final;

    if (current_covariance_.has_value()) {
        std::cout << "  - Covariance detected. Using STMPropagator for rigorous propagation.\n";
        
        // Setup initial cartesian state in InputFrame
        auto cart_initial = astdyn::propagation::keplerian_to_cartesian(kep_initial);
        Eigen::Vector<double, 6> x0;
        x0.head<3>() = cart_initial.position;
        x0.tail<3>() = cart_initial.velocity;

        // Setup STM Propagator
        auto integrator = std::make_unique<astdyn::propagation::RKF78Integrator>(
            settings_.initial_step, settings_.tolerance);
        
        // Force function from underlying propagator
        astdyn::propagation::STMPropagator stm_prop(
            std::move(integrator),
            [this](double t, const Eigen::Vector<double, 6>& x) {
                return propagator_->compute_derivatives(t, x);
            }
        );

        auto res = stm_prop.propagate(x0, kep_initial.epoch_mjd_tdb, target_mjd_tdb);
        pos_final = res.state.template head<3>();
        vel_final = res.state.template tail<3>();

        // Propagate Covariance: P(t) = Phi * P(t0) * Phi^T
        cov_final = res.stm * (*current_covariance_) * res.stm.transpose();
    } else {
        // Propagazione standard senza STM
        auto kep_final = propagator_->propagate_keplerian(kep_initial, target_mjd_tdb);
        auto cart_final = astdyn::propagation::keplerian_to_cartesian(kep_final);
        pos_final = cart_final.position;
        vel_final = cart_final.velocity;
    }
    
    std::cout << "  - Propagation finished." << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Statistiche
    std::ostringstream oss;
    oss << "Propagazione: " << current_epoch_mjd_ << " → " << target_mjd_tdb 
        << " MJD (" << (target_mjd_tdb - current_epoch_mjd_) << " giorni)\n"
        << "Tempo: " << duration.count() << " ms";
    if (cov_final) oss << " (STM Covariance propagated)";
    last_stats_ = oss.str();
    
    // Aggiorna elementi correnti per chiamate successive (getKeplerianElements)
    // Nota: se abbiamo usato STM, dovremmo riconvertire pos/vel in Kepleriani
    astdyn::propagation::CartesianElements cart_final_elem;
    cart_final_elem.position = pos_final;
    cart_final_elem.velocity = vel_final;
    cart_final_elem.epoch_mjd_tdb = target_mjd_tdb;
    cart_final_elem.gravitational_parameter = astdyn::constants::GMS;

    auto kep_final_updated = astdyn::propagation::cartesian_to_keplerian(cart_final_elem);
    current_elements_.semi_major_axis = kep_final_updated.semi_major_axis;
    current_elements_.eccentricity = kep_final_updated.eccentricity;
    current_elements_.inclination = kep_final_updated.inclination;
    current_elements_.longitude_asc_node = kep_final_updated.longitude_ascending_node;
    current_elements_.argument_perihelion = kep_final_updated.argument_perihelion;
    current_elements_.mean_anomaly = kep_final_updated.mean_anomaly;
    current_elements_.epoch_mjd_tdb = target_mjd_tdb;
    current_epoch_mjd_ = target_mjd_tdb;

    // Converti in ICRF se necessario
    if (current_frame_ == astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC) {
        auto result = eclipticToICRF(cart_final_elem, target_mjd_tdb);
        // Ruota anche la covarianza se presente
        if (cov_final) {
            // Matrice di rotazione Eclittica -> ICRF (3x3)
            const double cos_eps = std::cos(OrbitalConversions::OBLIQUITY_J2000);
            const double sin_eps = std::sin(OrbitalConversions::OBLIQUITY_J2000);
            Eigen::Matrix3d R;
            R << 1, 0, 0,
                 0, cos_eps, -sin_eps,
                 0, sin_eps, cos_eps;
            
            // Matrice di rotazione 6x6 (blocco diagonale)
            Eigen::Matrix<double, 6, 6> R6;
            R6.setZero();
            R6.template block<3, 3>(0, 0) = R;
            R6.template block<3, 3>(3, 3) = R;

            result.covariance = R6 * (*cov_final) * R6.transpose();
        }
        return result;
    } else {
        // Già in ICRF
        CartesianStateICRF result;
        result.position = pos_final;
        result.velocity = vel_final;
        result.epoch_mjd_tdb = target_mjd_tdb;
        result.covariance = cov_final;
        return result;
    }
}

FitResult AstDynWrapper::runFit(const std::string& rwo_filepath, bool verbose) {
    FitResult result;
    try {
        // Troviamo il file .eq1 originale (assumiamo sia nella stessa dir dell'oggetto se caricato)
        // Per ora richiediamo che loadFromEQ1 sia stato chiamato o passiamo un dummy
        // Implementazione via OrbitFitAPI
        std::string dummy_eq1 = ""; // Se non disponibile, l'API potrebbe fallire o richiedere input.
        
        // Poiché l'wrapper non salva il path dell'eq1, assumiamo che runFit
        // sia inteso per fittare l'orbita *già caricata*.
        // Ma l'API OrbitFitAPI::run_fit prende il path.
        
        result.message = "Fitting non ancora implementato via wrapper (usare OrbitFitAPI direttamente)";
        result.success = false;
        
        // Nota per l'utente: L'integrazione di runFit richiede il path del file .eq1 originale.
    } catch (const std::exception& e) {
        result.success = false;
        result.message = e.what();
    }
    return result;
}

astdyn::propagation::HighPrecisionPropagator::ObservationResult 
AstDynWrapper::calculateObservation(double target_mjd_tdb) {
    if (!initialized_) {
        std::cerr << "[AstDynWrapper:" << this << "] calculateObservation: ERROR - initialized_ is " << (initialized_ ? "TRUE" : "FALSE") << "!" << std::endl;
        throw std::runtime_error("Wrapper non inizializzato");
    }

    // Configura HighPrecisionPropagator con DE441 part-2 (lo standard validato)
    if (!high_prop_) {
        astdyn::propagation::HighPrecisionPropagator::Config prop_cfg;
        prop_cfg.de441_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
        prop_cfg.tolerance = 1e-13;
        high_prop_ = std::make_unique<astdyn::propagation::HighPrecisionPropagator>(prop_cfg);
        high_prop_->setPlanetaryEphemeris(ephemeris_);
    }

    // 1. Prepare elements
    astdyn::propagation::KeplerianElements kep;
    kep.epoch_mjd_tdb = current_elements_.epoch_mjd_tdb;
    kep.semi_major_axis = current_elements_.semi_major_axis;
    kep.eccentricity = current_elements_.eccentricity;
    kep.inclination = current_elements_.inclination;
    kep.longitude_ascending_node = current_elements_.longitude_asc_node;
    kep.argument_perihelion = current_elements_.argument_perihelion;
    kep.mean_anomaly = current_elements_.mean_anomaly;
    
    // Ensure gravitational parameter is set
    if (std::abs(current_elements_.gravitational_parameter) < 1e-18) {
        kep.gravitational_parameter = astdyn::constants::GMS;
    } else {
        kep.gravitational_parameter = current_elements_.gravitational_parameter;
    }

    double target_jd = target_mjd_tdb + 2400000.5;
    
    // IMPORTANTE: Poiché loadFromEQ1 applica convert_mean_to_osculating, 
    // gli elementi salvati sono in Equatorial/ICRF. Se usiamo setKeplerianElements,
    // usiamo il frame specificato.
    return high_prop_->calculateGeocentricObservation(
        kep, target_jd, current_frame_);
}

CartesianStateICRF AstDynWrapper::eclipticToICRF(
    const astdyn::propagation::CartesianElements& ecl_state,
    double epoch_mjd_tdb)
{
    // Obliquità eclittica J2000
    const double cos_eps = std::cos(OrbitalConversions::OBLIQUITY_J2000);
    const double sin_eps = std::sin(OrbitalConversions::OBLIQUITY_J2000);
    
    // Rotazione posizione
    const double x_ecl = ecl_state.position.x();
    const double y_ecl = ecl_state.position.y();
    const double z_ecl = ecl_state.position.z();
    
    Eigen::Vector3d pos_icrf;
    pos_icrf.x() = x_ecl;
    pos_icrf.y() = y_ecl * cos_eps - z_ecl * sin_eps;
    pos_icrf.z() = y_ecl * sin_eps + z_ecl * cos_eps;
    
    // Rotazione velocità
    const double vx_ecl = ecl_state.velocity.x();
    const double vy_ecl = ecl_state.velocity.y();
    const double vz_ecl = ecl_state.velocity.z();
    
    Eigen::Vector3d vel_icrf;
    vel_icrf.x() = vx_ecl;
    vel_icrf.y() = vy_ecl * cos_eps - vz_ecl * sin_eps;
    vel_icrf.z() = vy_ecl * sin_eps + vz_ecl * cos_eps;
    
    // Crea risultato
    CartesianStateICRF result;
    result.position = pos_icrf;
    result.velocity = vel_icrf;
    result.epoch_mjd_tdb = epoch_mjd_tdb;
    
    return result;
}

void AstDynWrapper::setSPKReader(std::shared_ptr<void> reader) {
    // No-op: AstDyn uses internal ephemeris
}

AstDynWrapper::ApparentState AstDynWrapper::getApparentStateGeocentric(double mjd_tdb) {
    auto obs = calculateObservation(mjd_tdb);
    ApparentState st;
    st.ra_deg = obs.ra_deg;
    st.dec_deg = obs.dec_deg;
    st.distance_au = obs.distance_au;
    st.position = obs.geocentric_position;
    return st;
}

astdyn::io::IOrbitParser::OrbitalElements AstDynWrapper::getKeplerianElements() const {
    return current_elements_;
}

std::string AstDynWrapper::getObjectName() const {
    return object_name_;
}

std::string AstDynWrapper::extractObjectNameFromEQ1(const std::string& filepath) {
    // 1. Prova a leggere dal file
    std::ifstream f(filepath);
    if (f.is_open()) {
        std::string line;
        while (std::getline(f, line)) {
            if (line.find("NAME") != std::string::npos && line.find("=") != std::string::npos) {
                 auto pos = line.find("=");
                 std::string name = line.substr(pos + 1);
                 // Trim
                 name.erase(0, name.find_first_not_of(" \t\r\n'\""));
                 name.erase(name.find_last_not_of(" \t\r\n'\"") + 1);
                 if (!name.empty()) return name;
            }
        }
    }
    
    // 2. Fallback: estrai dal nome del file (basename senza estensione)
    size_t last_slash = filepath.find_last_of("/\\");
    std::string filename = (last_slash == std::string::npos) ? filepath : filepath.substr(last_slash + 1);
    size_t last_dot = filename.find_last_of(".");
    if (last_dot != std::string::npos) {
        return filename.substr(0, last_dot);
    }
    return filename;
}

} // namespace ioccultcalc
