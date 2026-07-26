/**
 * @file test_external_validation.cpp
 * @brief Validazione contro riferimenti ESTERNI e indipendenti.
 *
 * Motivazione: i bug trovati il 2026-07-22 (rotazione terrestre invertita,
 * database osservatori vuoto con fallback silenzioso al geocentro, passo
 * dell'integratore che collassa) erano tutti SILENZIOSI e nessuno dei 170 test
 * esistenti li ha intercettati. Sono emersi solo confrontando i risultati con
 * riferimenti esterni: AstDyS per i residui, JPL Horizons per le posizioni,
 * i primi principi per la rotazione terrestre.
 *
 * Questi test rendono permanente quel confronto.
 *
 * I test che richiedono effemeridi o fixture vengono SALTATI (non falliti) se i
 * dati non sono disponibili, cosi' la suite resta eseguibile ovunque.
 */
#include <gtest/gtest.h>
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <astdyn/observations/ObservatoryDatabase.hpp>
#include <astdyn/observations/RWOReader.hpp>
#include <astdyn/orbit_determination/Residuals.hpp>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace astdyn;

namespace {

constexpr double AU_KM = 149597870.700;
constexpr double EPOCA_BK290 = 61200.0;

const char* ephemeris_path() { return std::getenv("ASTDYN_EPHEMERIS_PATH"); }

bool file_esiste(const std::string& p) { return std::ifstream(p).good(); }

physics::KeplerianStateTyped<core::ECLIPJ2000> orbita_bk290() {
    // Elementi AstDyS di 820987 (2015 BK290), epoca MJD 61200.0 TDT
    return physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(EPOCA_BK290),
        2.68777076554469, 0.103208119159221,
        11.8524777541786, 253.159103985645, 98.1398343280325, 333.016364934386);
}

AstDynEngine engine_standard() {
    AstDynEngine e;
    AstDynConfig c = e.config();
    c.verbose = false;
    c.integrator_type = IntegratorType::RKF78;
    c.tolerance = 1e-12;
    c.propagator_settings.include_planets = true;
    if (ephemeris_path()) c.ephemeris_file = ephemeris_path();
    c.ephemeris_type = EphemerisType::DE441;
    e.set_config(c);
    e.set_initial_orbit(orbita_bk290());
    return e;
}

} // namespace

// ===========================================================================
// 1. Rotazione terrestre — verifica dai PRIMI PRINCIPI
//    Non richiede dati esterni. Avrebbe intercettato il bug del segno di ERA.
// ===========================================================================

TEST(ValidazioneEsterna, RotazioneTerrestreGreenwich) {
    // Un punto sull'equatore a longitudine 0 (Greenwich) sta in ITRF in (R,0,0).
    // Dopo una rotazione terrestre di ERA, in GCRF deve trovarsi in
    //   (R cos ERA, +R sin ERA, 0)
    // Il segno di y e' il discriminante: se negativo, la Terra ruota al contrario.
    const double R = 6378137.0;
    Eigen::Vector3d greenwich_itrf(R, 0.0, 0.0);

    for (double mjd : {51544.5, 51544.75, 51545.0, 60000.0, 61200.0}) {
        auto t = time::EpochUTC::from_mjd(mjd);
        double era = coordinates::ReferenceFrame::computeERA(t, 0.0);
        auto M = coordinates::ReferenceFrame::itrf_to_j2000_simple(t);
        Eigen::Vector3d g = M * greenwich_itrf;

        EXPECT_NEAR(g.x(), R * std::cos(era), 1.0) << "MJD " << mjd;
        EXPECT_NEAR(g.y(), R * std::sin(era), 1.0) << "MJD " << mjd
            << " — segno invertito significa Terra che ruota al contrario";
        EXPECT_NEAR(g.z(), 0.0, 1.0) << "MJD " << mjd;
        EXPECT_NEAR(g.norm(), R, 1.0) << "la rotazione deve conservare il modulo";
    }
}

TEST(ValidazioneEsterna, ERAValoreNotoAJ2000) {
    // Valore IERS di riferimento: ERA = 280.46061838 gradi a J2000.0
    auto t = time::EpochUTC::from_mjd(51544.5);
    double era_deg = coordinates::ReferenceFrame::computeERA(t, 0.0) * 180.0 / M_PI;
    EXPECT_NEAR(era_deg, 280.46061838, 1e-4);
}

// ===========================================================================
// 2. Posizioni contro JPL Horizons
//    Oracoli a MJD INTERI esatti (un disallineamento di mezzo giorno
//    introdurrebbe ~5e-3 AU di errore spurio).
// ===========================================================================

TEST(ValidazioneEsterna, PosizioniBK290ControJPL) {
    if (!ephemeris_path()) GTEST_SKIP() << "ASTDYN_EPHEMERIS_PATH non impostata";

    struct Oracolo { double mjd; double x, y, z; };
    // JPL Horizons, 820987, eliocentrico eclittico J2000, AU
    const std::vector<Oracolo> oracoli = {
        {61200.0,  1.781080307462242, -1.615391586772313,  0.455977411927073},
        {60300.0, -1.848937426673913,  2.158158630411470, -0.502632346804024},
        {59000.0, -2.865115458643380, -0.453250895536448, -0.548392991206623},
        {57388.0, -2.867876932046533, -0.428057998980700, -0.550367845691710},
    };

    for (const auto& o : oracoli) {
        auto e = engine_standard();
        auto k = e.propagate_to(time::EpochTDB::from_mjd(o.mjd));
        auto p = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(k)
                     .position.to_eigen_si() / (AU_KM * 1000.0);
        double err = (p - Eigen::Vector3d(o.x, o.y, o.z)).norm();
        // 1e-5 AU = 1500 km ~ 1 arcsec a 2 AU: soglia larga, coglie i guasti veri
        EXPECT_LT(err, 1e-5) << "MJD " << o.mjd << " — scarto da JPL " << err << " AU";
    }
}

// ===========================================================================
// 3. Residui contro AstDyS — il test piu' importante
//    Copre l'intera catena: effemeridi, frame, tempo-luce, parallasse
//    topocentrica, propagazione. E' il confronto che ha trovato tutti i bug.
// ===========================================================================

TEST(ValidazioneEsterna, ResiduiBK290ControAstDyS) {
    if (!ephemeris_path()) GTEST_SKIP() << "ASTDYN_EPHEMERIS_PATH non impostata";

    std::string rwo;
    if (const char* d = std::getenv("ASTDYN_TEST_DATA")) {
        std::string cand = std::string(d) + "/820987.rwo";
        if (file_esiste(cand)) rwo = cand;
    }
    if (rwo.empty()) {
        for (const char* c : {"tests/data/820987.rwo", "../tests/data/820987.rwo",
                              "../../tests/data/820987.rwo"}) {
            if (file_esiste(c)) { rwo = c; break; }
        }
    }
    if (!file_esiste(rwo)) GTEST_SKIP() << "fixture mancante: " << rwo;

    // Il catalogo osservatori e' NECESSARIO: senza, la riduzione avviene dal
    // geocentro e i residui salgono da 0.199 a 1.35 arcsec. Il test si salta
    // invece di fallire, perche' un dato mancante non e' un errore del codice.
    // Ricerca nei percorsi convenzionali, come fa ioccultcalc.
    {
        std::string oc;
        if (const char* e = std::getenv("ASTDYN_OBSCODES")) {
            if (file_esiste(e)) oc = e;
        }
        if (oc.empty()) {
            const char* home = std::getenv("HOME");
            for (const std::string& c : {
                     std::string(home ? home : "") + "/.ioccultcalc/observatories/ObsCodes.txt",
                     std::string(home ? home : "") + "/.ioccultcalc/ObsCodes.txt"}) {
                if (file_esiste(c)) { oc = c; break; }
            }
        }
        if (oc.empty()) {
            GTEST_SKIP() << "catalogo osservatori non trovato: senza, la riduzione "
                            "avverrebbe dal geocentro e il confronto non avrebbe senso. "
                            "Installarlo con tools/ioccultcalc_setup.py --obscodes";
        }
        observations::ObservatoryDatabase::getInstance().loadFromMPCFile(oc);
    }

    auto e = engine_standard();
    auto cart0 = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(orbita_bk290());
    auto obs = observations::RWOReader::readFile(rwo);
    ASSERT_GT(obs.size(), 50u) << "il .rwo di riferimento deve contenere ~90 osservazioni";

    orbit_determination::ResidualCalculator<core::ECLIPJ2000> rc(e.getEphemeris(),
                                                                 e.propagator());
    auto res = rc.compute_residuals(obs, cart0);
    ASSERT_FALSE(res.empty());

    double somma = 0.0, somma_ra = 0.0, somma_dec = 0.0;
    for (const auto& r : res) {
        somma_ra  += r.residual_ra.to_arcsec();
        somma_dec += r.residual_dec.to_arcsec();
        somma     += std::abs(r.residual_ra.to_arcsec());
    }
    const double n = static_cast<double>(res.size());
    const double medio = somma / n;
    const double bias_ra = somma_ra / n, bias_dec = somma_dec / n;

    // Riferimento AstDyS sullo stesso insieme: |residuo| medio 0.169 arcsec.
    // Misura del 2026-07-22 dopo le correzioni: 0.199 arcsec.
    // Soglia a 0.5: coglie i guasti sistematici senza essere fragile.
    EXPECT_LT(medio, 0.5) << "residuo medio " << medio
        << " arcsec (riferimento AstDyS 0.169). Un valore molto piu' alto indica "
           "un errore sistematico nel modello di osservazione.";

    // I bias devono essere prossimi a zero: un offset costante e' la firma
    // caratteristica di un termine mancante o con segno errato.
    EXPECT_LT(std::abs(bias_ra),  0.3) << "bias sistematico in RA: "  << bias_ra;
    EXPECT_LT(std::abs(bias_dec), 0.3) << "bias sistematico in Dec: " << bias_dec;
}

// ===========================================================================
// 4. Database osservatori — il fallback al geocentro non deve essere silenzioso
// ===========================================================================

TEST(ValidazioneEsterna, DatabaseOsservatoriPopolato) {
    const char* obscodes = std::getenv("ASTDYN_OBSCODES");
    if (!obscodes || !file_esiste(obscodes)) {
        GTEST_SKIP() << "ASTDYN_OBSCODES non disponibile";
    }
    auto& db = observations::ObservatoryDatabase::getInstance();
    size_t n = db.loadFromMPCFile(obscodes);
    EXPECT_GT(n, 2000u) << "il catalogo MPC contiene oltre 2000 codici; "
                           "un numero piccolo significa che si sta usando la "
                           "lista ridotta e molte osservazioni cadranno sul geocentro";

    // codici molto usati nell'astrometria di asteroidi
    for (const char* c : {"691", "F51", "G96", "703", "W68", "D04"}) {
        EXPECT_TRUE(db.getObservatory(c).has_value())
            << "codice " << c << " assente: le sue osservazioni userebbero il "
               "geocentro, perdendo la parallasse topocentrica";
    }
}

// ===========================================================================
// 5. Integratori — concordanza reciproca
//    Metodi indipendenti devono concordare; una divergenza segnala un guasto.
// ===========================================================================

TEST(ValidazioneEsterna, IntegratoriConcordanti) {
    if (!ephemeris_path()) GTEST_SKIP() << "ASTDYN_EPHEMERIS_PATH non impostata";

    const double target = EPOCA_BK290 - 365.0;   // un anno indietro

    auto propaga = [&](IntegratorType tipo) {
        AstDynEngine e;
        AstDynConfig c = e.config();
        c.verbose = false;
        c.integrator_type = tipo;
        c.tolerance = 1e-12;
        c.propagator_settings.include_planets = true;
        c.ephemeris_file = ephemeris_path();
        c.ephemeris_type = EphemerisType::DE441;
        e.set_config(c);
        e.set_initial_orbit(orbita_bk290());
        auto k = e.propagate_to(time::EpochTDB::from_mjd(target));
        return propagation::keplerian_to_cartesian<core::ECLIPJ2000>(k)
                   .position.to_eigen_si() / (AU_KM * 1000.0);
    };

    auto rif = propaga(IntegratorType::RKF78);
    for (auto tipo : {IntegratorType::RK4, IntegratorType::GAUSS,
                      IntegratorType::GRKN64, IntegratorType::AAS}) {
        double err = (propaga(tipo) - rif).norm();
        EXPECT_LT(err, 1e-6) << "scarto da RKF78 su un anno: " << err << " AU";
    }
}

TEST(ValidazioneEsterna, SABA4RifiutatoEsplicitamente) {
    // SABA4 e' difettoso (errore 1.60 AU su 0.1 giorni): deve fallire in modo
    // esplicito, non restituire silenziosamente risultati errati.
    AstDynEngine e;
    AstDynConfig c = e.config();
    c.integrator_type = IntegratorType::SABA4;
    EXPECT_THROW(e.set_config(c), std::runtime_error);
}
