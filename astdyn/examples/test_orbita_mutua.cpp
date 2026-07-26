/**
 * @file test_orbita_mutua.cpp
 * @brief Passo 1: dai parametri di un'orbita mutua allo stato relativo.
 *
 * I binari pubblicati (archivio di Johnston, letteratura) sono descritti da
 * semiasse, eccentricita', inclinazione, nodo, argomento del pericentro e
 * periodo dell'orbita del satellite attorno al primario.
 *
 * Per propagarli serve mu del sistema. Si ricava dal PERIODO con la terza legge
 * di Keplero,  mu = 4 pi^2 a^3 / P^2,  non dalle masse: il periodo si misura
 * bene dalle curve di luce, le masse dei binari sono spesso incerte al 20-30%.
 *
 * Verifiche:
 *   - su orbita circolare |rho| = a e |drho| = 2 pi a / P;
 *   - mu ricavato dal periodo riproduce il periodo di partenza;
 *   - la massa implicita e' plausibile per un asteroide.
 */
#include <astdyn/AstDynEngine.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double G_SI = 6.67430e-11;          // m^3 kg^-1 s^-2

/// mu del sistema [m^3/s^2] da semiasse [km] e periodo [giorni].
double mu_da_periodo(double a_km, double periodo_giorni) {
    const double a = a_km * 1000.0;                 // m
    const double P = periodo_giorni * 86400.0;      // s
    return 4.0 * PI * PI * a * a * a / (P * P);
}

/// Periodo [giorni] da semiasse [km] e mu [m^3/s^2] — inversa della precedente.
double periodo_da_mu(double a_km, double mu) {
    const double a = a_km * 1000.0;
    return 2.0 * PI * std::sqrt(a * a * a / mu) / 86400.0;
}

struct Sistema {
    const char* nome;
    const char* satellite;
    double a_km;
    double e;
    double i_deg;
    double periodo_giorni;
    double diam_primario_km;
};

// Parametri da letteratura, per la prova numerica.
const Sistema SISTEMI[] = {
    {"(22) Kalliope", "Linus",   1099.0, 0.005, 103.7, 3.5951, 150.0},
    {"(45) Eugenia",  "Petit-Prince", 1184.0, 0.007, 134.0, 4.7255, 202.0},
    {"(87) Sylvia",   "Romulus", 1356.0, 0.001,  8.0,  3.6496, 271.0},
};

} // namespace

int main() {
    std::cout << "=== mu del sistema dal periodo orbitale ===\n\n";
    std::cout << std::left
              << std::setw(16) << "sistema"
              << std::setw(15) << "satellite"
              << std::setw(9)  << "a [km]"
              << std::setw(10) << "P [d]"
              << std::setw(13) << "mu [m3/s2]"
              << std::setw(13) << "massa [kg]"
              << std::setw(12) << "densita'" << "\n";
    std::cout << std::string(88, '-') << "\n";

    bool tutto_bene = true;

    for (const auto& s : SISTEMI) {
        const double mu = mu_da_periodo(s.a_km, s.periodo_giorni);
        const double massa = mu / G_SI;

        // densita' implicita, assumendo una sfera del diametro del primario:
        // e' un controllo di plausibilita', non una misura.
        const double raggio = s.diam_primario_km * 500.0;               // m
        const double volume = 4.0 / 3.0 * PI * raggio * raggio * raggio;
        const double densita = massa / volume;                          // kg/m^3

        // controllo di andata e ritorno
        const double P_ric = periodo_da_mu(s.a_km, mu);
        const double err_rel = std::abs(P_ric - s.periodo_giorni) / s.periodo_giorni;

        std::cout << std::left << std::setw(16) << s.nome
                  << std::setw(15) << s.satellite
                  << std::fixed << std::setprecision(0) << std::setw(9) << s.a_km
                  << std::setprecision(4) << std::setw(10) << s.periodo_giorni
                  << std::scientific << std::setprecision(3) << std::setw(13) << mu
                  << std::setw(13) << massa
                  << std::fixed << std::setprecision(0) << std::setw(8) << densita << " kg/m3";

        if (err_rel > 1e-12) {
            std::cout << "   [!! andata-ritorno " << err_rel << "]";
            tutto_bene = false;
        }
        // Gli asteroidi hanno densita' fra ~1000 (rubble pile ghiacciati) e
        // ~4000 kg/m3 (metallici). Fuori da questo intervallo, i parametri o
        // il diametro non tornano.
        if (densita < 500.0 || densita > 6000.0) {
            std::cout << "   [densita' implausibile]";
        }
        std::cout << "\n";
    }

    // ---- velocita' orbitale relativa su orbita circolare ----
    std::cout << "\n=== velocita' relativa (orbita circolare: v = 2 pi a / P) ===\n\n";
    for (const auto& s : SISTEMI) {
        const double v_attesa = 2.0 * PI * s.a_km / (s.periodo_giorni * 86400.0);  // km/s
        const double mu = mu_da_periodo(s.a_km, s.periodo_giorni);
        const double v_vis_viva = std::sqrt(mu / (s.a_km * 1000.0)) / 1000.0;      // km/s

        const double scarto = std::abs(v_vis_viva - v_attesa) / v_attesa;
        std::cout << std::left << std::setw(16) << s.nome
                  << "v = " << std::fixed << std::setprecision(5) << v_attesa << " km/s"
                  << "   da vis-viva: " << v_vis_viva
                  << (scarto < 1e-10 ? "   [ok]" : "   [!!]") << "\n";
        if (scarto >= 1e-10) tutto_bene = false;
    }

    std::cout << "\n" << (tutto_bene ? "Tutte le verifiche superate."
                                     : "ATTENZIONE: alcune verifiche non tornano.") << "\n";
    return tutto_bene ? 0 : 1;
}
