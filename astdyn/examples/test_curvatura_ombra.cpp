/**
 * @file test_curvatura_ombra.cpp
 * @brief Quanto costa approssimare il moto dell'ombra con una retta?
 *
 * `OccultationMapper::calculate_geopoint_at_epoch` calcola la posizione
 * dell'asse dell'ombra come
 *      xi(t) = xi_ca + dxi_dt * t
 * cioe' per estrapolazione lineare dalle derivate al massimo avvicinamento.
 * Il moto vero e' curvo: il formato Occult4 prevede infatti termini fino al
 * terzo ordine, che nelle nostre esportazioni sono a zero.
 *
 * Questo programma confronta, per un evento reale, la posizione estrapolata
 * linearmente con quella ricalcolata dalla geometria a istanti crescenti dal
 * TCA, e riporta lo scarto in chilometri sul piano fondamentale.
 *
 * La domanda cui risponde: l'approssimazione vale pochi chilometri, come dice
 * il commento nel codice, o decine?
 */
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/astrometry/ChebyshevEphemerisManager.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;

namespace {

constexpr double DEG = 3.14159265358979323846 / 180.0;
constexpr double AU_M = 1.495978707e11;

/// Coordinate del piano fondamentale: proiezione della congiungente
/// Terra-asteroide sul piano perpendicolare alla direzione della stella.
/// Replica di compute_fundamental_plane_geometry, che e' privata.
void piano_fondamentale(double star_ra_deg, double star_dec_deg,
                        double ast_ra_deg, double ast_dec_deg, double dist_au,
                        double& xi_m, double& eta_m) {
    const double as = star_ra_deg * DEG, ds = star_dec_deg * DEG;
    const double a  = ast_ra_deg  * DEG, d  = ast_dec_deg  * DEG;
    const double r  = dist_au * AU_M;

    // versore verso l'asteroide
    const double ux = std::cos(d) * std::cos(a);
    const double uy = std::cos(d) * std::sin(a);
    const double uz = std::sin(d);

    // versore verso la stella e assi est/nord del piano fondamentale
    const double sx = std::cos(ds) * std::cos(as);
    const double sy = std::cos(ds) * std::sin(as);
    const double sz = std::sin(ds);

    const double ex = -std::sin(as),          ey = std::cos(as),           ez = 0.0;
    const double nx = -std::sin(ds)*std::cos(as),
                 ny = -std::sin(ds)*std::sin(as),
                 nz =  std::cos(ds);

    // componente della posizione dell'asteroide perpendicolare alla stella
    const double px = r*ux, py = r*uy, pz = r*uz;
    const double lungo = px*sx + py*sy + pz*sz;
    const double qx = px - lungo*sx, qy = py - lungo*sy, qz = pz - lungo*sz;

    xi_m  = qx*ex + qy*ey + qz*ez;
    eta_m = qx*nx + qy*ny + qz*nz;
}

} // namespace

int main(int argc, char** argv) {
    const char* eph = std::getenv("ASTDYN_EPHEMERIS_PATH");
    if (!eph) { std::cout << "serve ASTDYN_EPHEMERIS_PATH\n"; return 0; }

    // Evento 820987 del 2026-07-27, TCA a MJD 61248.0472 (01:07:59 UT)
    const double mjd_tca = 61248.0 + 1.1325791 / 24.0;
    const double star_ra_deg  = 22.84268511 * 15.0;   // il file XML ha ore
    const double star_dec_deg = 11.60197462;

    AstDynConfig cfg;
    cfg.verbose = false;
    cfg.integrator_type = IntegratorType::RKF78;
    cfg.tolerance = 1e-11;
    cfg.propagator_settings.include_planets = true;
    cfg.ephemeris_file = eph;
    cfg.ephemeris_type = EphemerisType::DE441;

    auto orbita = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(61200.0),
        2.68777076554469, 0.103208119159221, 11.8524777541786,
        253.159103985645, 98.1398343280325, 333.016364934386);

    astrometry::ChebyshevEphemerisManager manager(cfg);
    manager.add_asteroid("820987", orbita,
                         time::EpochTDB::from_mjd(mjd_tca - 0.2),
                         time::EpochTDB::from_mjd(mjd_tca + 0.2));

    // geometria al TCA e derivate prime
    auto valuta = [&](double mjd, double& xi, double& eta, double& ra, double& dec) {
        const auto t = time::EpochTDB::from_mjd(mjd);
        const auto [pos, vel] = manager.evaluate_full("820987", t);
        const auto [r, d, dist] = pos;
        ra = r; dec = d;
        piano_fondamentale(star_ra_deg, star_dec_deg, r, d, dist, xi, eta);
    };

    double xi0, eta0, ra0, dec0;
    valuta(mjd_tca, xi0, eta0, ra0, dec0);

    // derivate per differenze centrate su un minuto
    const double h = 1.0 / 1440.0;
    double xip, etap, xim, etam, rr, dd;
    valuta(mjd_tca + h, xip, etap, rr, dd);
    valuta(mjd_tca - h, xim, etam, rr, dd);
    const double dxi  = (xip - xim)  / (2.0 * h * 86400.0);   // m/s
    const double deta = (etap - etam) / (2.0 * h * 86400.0);

    std::cout << "=== moto dell'ombra: retta contro geometria vera ===\n\n";
    std::cout << "TCA          MJD " << std::fixed << std::setprecision(6) << mjd_tca << "\n";
    std::cout << "xi, eta      " << std::setprecision(1) << xi0/1000.0 << ", "
              << eta0/1000.0 << " km\n";
    std::cout << "velocita'    " << std::setprecision(3)
              << std::sqrt(dxi*dxi + deta*deta)/1000.0 << " km/s\n\n";

    std::cout << std::left << std::setw(12) << "dt [min]"
              << std::setw(16) << "scarto xi [km]"
              << std::setw(16) << "scarto eta [km]"
              << std::setw(14) << "totale [km]" << "\n";
    std::cout << std::string(58, '-') << "\n";

    double peggiore = 0.0;
    for (double minuti : {-120.0, -90.0, -60.0, -30.0, -15.0, 15.0, 30.0, 60.0, 90.0, 120.0}) {
        const double dt_s = minuti * 60.0;
        const double mjd = mjd_tca + minuti / 1440.0;

        double xi_vero, eta_vero, r, d;
        valuta(mjd, xi_vero, eta_vero, r, d);

        const double xi_lin  = xi0  + dxi  * dt_s;
        const double eta_lin = eta0 + deta * dt_s;

        const double sx = (xi_vero - xi_lin) / 1000.0;
        const double sy = (eta_vero - eta_lin) / 1000.0;
        const double tot = std::sqrt(sx*sx + sy*sy);
        peggiore = std::max(peggiore, tot);

        std::cout << std::left << std::setw(12) << std::setprecision(0) << minuti
                  << std::setw(16) << std::setprecision(2) << sx
                  << std::setw(16) << sy
                  << std::setw(14) << tot << "\n";
    }

    std::cout << "\nscarto massimo su +/-2 ore: " << std::setprecision(1)
              << peggiore << " km\n";
    std::cout << "(la larghezza dell'ombra di questo corpo e' 1.25 km)\n\n";

    if (peggiore < 5.0)
        std::cout << "L'approssimazione lineare regge: pochi km, come dice il commento.\n";
    else
        std::cout << "L'approssimazione lineare NON regge: lo scarto e' molte volte\n"
                  << "la larghezza dell'ombra, e sposta i siti di osservazione.\n";
    return 0;
}
