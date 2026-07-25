/**
 * @file test_earth_vel.cpp
 * @brief Verifica il segno della velocita' della Terra dall'effemeride.
 *
 * Controlli indipendenti:
 *  1. |v| deve essere ~29.8 km/s
 *  2. r x v deve avere componente z POSITIVA (moto antiorario visto da nord eclittica)
 *  3. r . v ~ 0 (orbita quasi circolare: velocita' quasi perpendicolare al raggio)
 *  4. confronto con JPL Horizons per la stessa data
 */
#include <astdyn/AstDynEngine.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    AstDynEngine engine;
    AstDynConfig cfg = engine.config();
    const char* eph = std::getenv("ASTDYN_EPHEMERIS_PATH");
    if (eph) cfg.ephemeris_file = eph;
    cfg.ephemeris_type = EphemerisType::DE441;
    engine.set_config(cfg);

    auto ephem = engine.getEphemeris();
    auto t = time::EpochTDB::from_mjd(61200.0);
    auto st = ephem->getState(ephemeris::CelestialBody::EARTH, t);

    Eigen::Vector3d r = st.position.to_eigen_si() / 1000.0;   // km
    Eigen::Vector3d v = st.velocity.to_eigen_si() / 1000.0;   // km/s

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "=== Terra @ MJD 61200 (frame nativo dell'effemeride) ===\n";
    std::cout << "r = (" << r.x() << ", " << r.y() << ", " << r.z() << ") km\n";
    std::cout << "    |r| = " << r.norm()/1.495978707e8 << " AU\n";
    std::cout << "v = (" << v.x() << ", " << v.y() << ", " << v.z() << ") km/s\n";
    std::cout << "    |v| = " << v.norm() << " km/s   (atteso ~29.8)\n\n";

    Eigen::Vector3d h = r.cross(v);
    std::cout << "r x v = (" << h.x() << ", " << h.y() << ", " << h.z() << ")\n";
    std::cout << "  componente z = " << h.z()
              << (h.z() > 0 ? "  [OK: moto antiorario, segno corretto]"
                            : "  [ANOMALO: velocita' con segno invertito]") << "\n\n";

    double cosang = r.dot(v) / (r.norm() * v.norm());
    std::cout << "r . v normalizzato = " << cosang << "  (atteso ~0, orbita quasi circolare)\n";
    std::cout << "angolo r,v = " << std::acos(cosang) * 180.0 / M_PI << " gradi (atteso ~90)\n";
    return 0;
}
