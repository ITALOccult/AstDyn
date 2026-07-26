/**
 * @file test_satellite.cpp
 * @brief Verifica del vettore relativo di un satellite sull'orbita mutua.
 *
 * Controlli, dal piu' semplice al piu' stringente:
 *   1. orbita circolare: |rho| costante e pari al semiasse;
 *   2. periodo: dopo un periodo esatto la posizione si ripete;
 *   3. velocita': la derivata numerica coincide con 2 pi a / P;
 *   4. eccentrica: pericentro a(1-e) e apocentro a(1+e);
 *   5. inclinazione: l'ampiezza in z corrisponde a  a sin(i);
 *   6. piano di riferimento: eclittico ed equatoriale differiscono di ~23.4 gradi.
 */
#include <astdyn/astrometry/OrbitingChebyshevEphemeris.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;
using astrometry::OrbitingChebyshevEphemeris;
using Orbita = OrbitingChebyshevEphemeris::OrbitaMutua;
using Piano = OrbitingChebyshevEphemeris::PianoRiferimento;

namespace {

constexpr double PI = 3.14159265358979323846;
bool tutto_bene = true;

void verifica(const char* cosa, double ottenuto, double atteso, double tolleranza,
              const char* unita = "") {
    const double scarto = std::abs(ottenuto - atteso);
    const bool ok = scarto <= tolleranza;
    if (!ok) tutto_bene = false;
    std::cout << "  " << std::left << std::setw(38) << cosa
              << std::fixed << std::setprecision(6)
              << std::setw(14) << ottenuto << " atteso " << std::setw(14) << atteso
              << unita << (ok ? "  [ok]" : "  [!!]") << "\n";
}

} // namespace

int main() {
    const double MJD0 = 61200.0;

    // ---- 1, 2, 3: orbita circolare equatoriale nel piano eclittico ----
    std::cout << "=== orbita circolare (a=1099 km, P=3.5951 d, i=0) ===\n";
    Orbita circ;
    circ.a_km = 1099.0;
    circ.e = 0.0;
    circ.i_deg = 0.0;
    circ.node_deg = 0.0;
    circ.peri_deg = 0.0;
    circ.M_deg = 0.0;
    circ.period_days = 3.5951;
    circ.epoch = time::EpochTDB::from_mjd(MJD0);
    circ.piano = Piano::Eclittico;

    // il modulo deve restare pari al semiasse a ogni istante
    double max_scarto_raggio = 0.0;
    for (int k = 0; k <= 20; ++k) {
        const double t = MJD0 + circ.period_days * k / 20.0;
        const auto r = OrbitingChebyshevEphemeris::vettore_relativo(
            circ, time::EpochTDB::from_mjd(t));
        max_scarto_raggio = std::max(max_scarto_raggio, std::abs(r.norm() - circ.a_km));
    }
    verifica("modulo costante (scarto max)", max_scarto_raggio, 0.0, 1e-6, " km");

    // dopo un periodo esatto la posizione si ripete
    const auto r0 = OrbitingChebyshevEphemeris::vettore_relativo(
        circ, time::EpochTDB::from_mjd(MJD0));
    const auto rP = OrbitingChebyshevEphemeris::vettore_relativo(
        circ, time::EpochTDB::from_mjd(MJD0 + circ.period_days));
    verifica("ritorno dopo un periodo", (rP - r0).norm(), 0.0, 1e-6, " km");

    // all'epoca, con M=peri=node=0, il satellite sta sull'asse x a distanza a
    verifica("posizione iniziale x", r0.x(), circ.a_km, 1e-9, " km");
    verifica("posizione iniziale y", r0.y(), 0.0, 1e-9, " km");

    // velocita' per differenze finite
    const double dt_giorni = 1.0e-4;
    const auto rp = OrbitingChebyshevEphemeris::vettore_relativo(
        circ, time::EpochTDB::from_mjd(MJD0 + dt_giorni));
    const auto rm = OrbitingChebyshevEphemeris::vettore_relativo(
        circ, time::EpochTDB::from_mjd(MJD0 - dt_giorni));
    const double v = (rp - rm).norm() / (2.0 * dt_giorni * 86400.0);   // km/s
    const double v_attesa = 2.0 * PI * circ.a_km / (circ.period_days * 86400.0);
    verifica("velocita' orbitale", v, v_attesa, 1e-6, " km/s");

    // ---- 4: orbita eccentrica ----
    std::cout << "\n=== orbita eccentrica (e=0.3) ===\n";
    Orbita ecc = circ;
    ecc.e = 0.3;
    const auto r_peri = OrbitingChebyshevEphemeris::vettore_relativo(
        ecc, time::EpochTDB::from_mjd(MJD0));                             // M=0 -> pericentro
    const auto r_apo = OrbitingChebyshevEphemeris::vettore_relativo(
        ecc, time::EpochTDB::from_mjd(MJD0 + ecc.period_days / 2.0));      // M=pi -> apocentro
    verifica("distanza al pericentro", r_peri.norm(), ecc.a_km * (1.0 - ecc.e), 1e-6, " km");
    verifica("distanza all'apocentro", r_apo.norm(), ecc.a_km * (1.0 + ecc.e), 1e-6, " km");

    // ---- 5: inclinazione ----
    std::cout << "\n=== inclinazione (i=30 gradi) ===\n";
    Orbita inc = circ;
    inc.i_deg = 30.0;
    double z_max = 0.0;
    for (int k = 0; k <= 200; ++k) {
        const double t = MJD0 + inc.period_days * k / 200.0;
        const auto r = OrbitingChebyshevEphemeris::vettore_relativo(
            inc, time::EpochTDB::from_mjd(t));
        z_max = std::max(z_max, std::abs(r.z()));
    }
    verifica("escursione massima in z", z_max, inc.a_km * std::sin(30.0 * PI / 180.0),
             0.5, " km");

    // ---- 6: piano di riferimento ----
    std::cout << "\n=== piano di riferimento ===\n";
    Orbita eq = circ;  eq.piano = Piano::Equatoriale;
    Orbita ec = circ;  ec.piano = Piano::Eclittico;
    const auto r_eq = OrbitingChebyshevEphemeris::vettore_relativo(
        eq, time::EpochTDB::from_mjd(MJD0 + 0.9));
    const auto r_ec = OrbitingChebyshevEphemeris::vettore_relativo(
        ec, time::EpochTDB::from_mjd(MJD0 + 0.9));
    const double angolo = std::acos(r_eq.dot(r_ec) / (r_eq.norm() * r_ec.norm()))
                          * 180.0 / PI;
    std::cout << "  angolo fra le due interpretazioni: "
              << std::fixed << std::setprecision(3) << angolo << " gradi\n";
    std::cout << "  (per un'orbita nel piano di riferimento coincide con l'obliquita', 23.44)\n";
    if (angolo < 0.001) {
        std::cout << "  [!!] nessuna differenza: la conversione di piano non viene applicata\n";
        tutto_bene = false;
    }

    // ---- mu del sistema ----
    std::cout << "\n=== parametro gravitazionale ===\n";
    const double mu = OrbitingChebyshevEphemeris::mu_dal_periodo(1099.0, 3.5951);
    const double massa = mu * 1e9 / 6.67430e-11;    // km^3/s^2 -> m^3/s^2 -> kg
    std::cout << "  mu    = " << std::scientific << std::setprecision(4) << mu << " km3/s2\n";
    std::cout << "  massa = " << massa << " kg   (Kalliope: ~8.1e18 atteso)\n";
    if (massa < 5e18 || massa > 1.2e19) {
        std::cout << "  [!!] massa fuori dall'intervallo atteso\n";
        tutto_bene = false;
    }

    std::cout << "\n" << (tutto_bene ? "Tutte le verifiche superate."
                                     : "ATTENZIONE: alcune verifiche non tornano.") << "\n";
    return tutto_bene ? 0 : 1;
}
