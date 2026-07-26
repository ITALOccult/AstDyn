/**
 * @file OrbitingChebyshevEphemeris.cpp
 * @brief Effemeride di un satellite propagato sulla sua orbita mutua.
 */

#include "astdyn/astrometry/OrbitingChebyshevEphemeris.hpp"
#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace astdyn::astrometry {

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double GRADI_A_RAD = PI / 180.0;
constexpr double SECONDI_AL_GIORNO = 86400.0;

/// Risolve l'equazione di Keplero  E - e sin E = M  con Newton-Raphson.
/// Ritorna l'anomalia eccentrica in radianti.
double risolvi_keplero(double M, double e) {
    // normalizzazione in [-pi, pi] per una partenza vicina alla soluzione
    M = std::fmod(M, 2.0 * PI);
    if (M >  PI) M -= 2.0 * PI;
    if (M < -PI) M += 2.0 * PI;

    // stima iniziale: per e piccola E ~ M; per e grande conviene partire da pi
    double E = (e < 0.8) ? M : PI;

    for (int i = 0; i < 50; ++i) {
        const double f  = E - e * std::sin(E) - M;
        const double fp = 1.0 - e * std::cos(E);
        const double dE = f / fp;
        E -= dE;
        if (std::abs(dE) < 1e-14) break;
    }
    return E;
}

} // namespace

double OrbitingChebyshevEphemeris::mu_dal_periodo(double a_km, double period_days) {
    if (a_km <= 0.0 || period_days <= 0.0) {
        throw std::invalid_argument(
            "OrbitingChebyshevEphemeris: semiasse e periodo dell'orbita mutua devono "
            "essere positivi (a=" + std::to_string(a_km) + " km, P=" +
            std::to_string(period_days) + " d)");
    }
    const double P = period_days * SECONDI_AL_GIORNO;      // s
    return 4.0 * PI * PI * a_km * a_km * a_km / (P * P);   // km^3/s^2
}

Eigen::Vector3d OrbitingChebyshevEphemeris::vettore_relativo(
    const OrbitaMutua& orb, time::EpochTDB t)
{
    const double mu = mu_dal_periodo(orb.a_km, orb.period_days);

    // moto medio e anomalia media all'istante richiesto
    const double n = std::sqrt(mu / (orb.a_km * orb.a_km * orb.a_km));   // rad/s
    const double dt = (t.mjd() - orb.epoch.mjd()) * SECONDI_AL_GIORNO;   // s
    const double M = orb.M_deg * GRADI_A_RAD + n * dt;

    const double E = risolvi_keplero(M, orb.e);

    // posizione nel piano orbitale
    const double cosE = std::cos(E), sinE = std::sin(E);
    const double fattore = std::sqrt(1.0 - orb.e * orb.e);
    const double x_p = orb.a_km * (cosE - orb.e);
    const double y_p = orb.a_km * fattore * sinE;

    // rotazione al frame di riferimento: R_z(-nodo) R_x(-i) R_z(-peri)
    const double om = orb.peri_deg * GRADI_A_RAD;
    const double On = orb.node_deg * GRADI_A_RAD;
    const double in = orb.i_deg    * GRADI_A_RAD;

    const double co = std::cos(om), so = std::sin(om);
    const double cO = std::cos(On), sO = std::sin(On);
    const double ci = std::cos(in), si = std::sin(in);

    Eigen::Vector3d r;
    r.x() = (cO * co - sO * so * ci) * x_p + (-cO * so - sO * co * ci) * y_p;
    r.y() = (sO * co + cO * so * ci) * x_p + (-sO * so + cO * co * ci) * y_p;
    r.z() = (so * si)                * x_p + ( co * si)                * y_p;

    // I parametri pubblicati sono di norma riferiti all'equatore J2000, mentre
    // la libreria lavora in eclittica: la conversione va fatta esplicitamente,
    // altrimenti l'orbita risulta ruotata di ~23.4 gradi e le tracce finiscono
    // altrove in modo credibile.
    if (orb.piano == PianoRiferimento::Equatoriale) {
        const auto R = coordinates::ReferenceFrame::get_rotation<core::GCRF, core::ECLIPJ2000>(t);
        r = R.to_eigen() * r;
    }
    return r;
}

OrbitingChebyshevEphemeris::OrbitingChebyshevEphemeris(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& primary_elements,
    const OrbitaMutua& orbita,
    time::EpochTDB start_time,
    time::EpochTDB end_time,
    const astdyn::AstDynConfig& config,
    int degree)
    : start_epoch_(start_time), end_epoch_(end_time)
{
    jd_start_ = start_time.jd();
    jd_end_   = end_time.jd();

    if (jd_end_ <= jd_start_) {
        throw std::invalid_argument(
            "OrbitingChebyshevEphemeris: end_time deve seguire start_time");
    }

    // Un motore dedicato per propagare l'orbita eliocentrica del primario.
    // Costruito una volta sola: ricrearlo per ogni campione costerebbe caro.
    auto motore = std::make_shared<AstDynEngine>();
    motore->set_config(config);
    motore->set_initial_orbit(primary_elements);

    // Fornitore di posizione: primario propagato piu' vettore relativo.
    auto posizione_satellite = [motore, orbita](double et_tdb) -> Eigen::Vector3d {
        const double jd = et_tdb / SECONDI_AL_GIORNO + constants::JD2000;
        const time::EpochTDB t = time::EpochTDB::from_jd(jd);

        auto kep = motore->propagate_to(t);
        auto cart = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(kep);
        const Eigen::Vector3d r_primario = cart.position.to_eigen_si() / 1000.0;  // km

        return r_primario + vettore_relativo(orbita, t);
    };

    // Segmenti giornalieri, come per gli altri corpi. Un satellite con periodo
    // di qualche giorno percorre una frazione modesta di orbita in un giorno,
    // quindi un polinomio di grado 12 la descrive ampiamente.
    double current_jd = jd_start_;
    while (current_jd < jd_end_) {
        const double next_jd = std::min(current_jd + 1.0, jd_end_);
        const double durata = next_jd - current_jd;
        const auto centro = time::EpochTDB::from_jd(current_jd + durata / 2.0);

        segments_.push_back(
            catalog::fit_chebyshev_da_posizione(posizione_satellite, centro, durata, degree));
        current_jd = next_jd;
    }
}

std::tuple<double, double, double>
OrbitingChebyshevEphemeris::evaluate(time::EpochTDB epoch) const {
    return get_segment(epoch).evaluate_all(epoch.jd());
}

std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>>
OrbitingChebyshevEphemeris::evaluate_full(time::EpochTDB epoch) const {
    return get_segment(epoch).evaluate_full(epoch.jd());
}

const catalog::ChebyshevSegment&
OrbitingChebyshevEphemeris::get_segment(time::EpochTDB epoch) const {
    const double jd = epoch.jd();
    if (jd < jd_start_ || jd > jd_end_) {
        throw std::out_of_range("OrbitingChebyshevEphemeris: epoca fuori intervallo");
    }
    size_t idx = static_cast<size_t>(std::floor(jd - jd_start_));
    if (idx >= segments_.size()) idx = segments_.size() - 1;
    return segments_[idx];
}

} // namespace astdyn::astrometry
