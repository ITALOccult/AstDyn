/**
 * @file test_nbody_frame.cpp
 * @brief Verifica che la rotazione GCRF->ecliptic usata dal N-body porti un
 *        corpo dell'eclittica (es. la Terra) a latitudine eclittica ~0.
 *        Non richiede effemeride: usa la proprieta' geometrica invariante.
 *        La validazione di AMPIEZZA end-to-end e' nella prova su 2015 BK290.
 */
#include <gtest/gtest.h>
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>

using namespace astdyn;

// Una posizione realmente nell'eclittica, espressa in GCRF (equatoriale):
// in eclittica ha lat 0; vista in GCRF acquista lat = obliquita'.
TEST(NBodyFrame, EclipticBodyMapsToZeroLatitude) {
    const auto t = time::EpochTDB::from_mjd(60000.0);
    const Matrix3d ecl2eq = coordinates::ReferenceFrame::ecliptic_to_j2000(t);
    const Matrix3d eq2ecl = coordinates::ReferenceFrame::j2000_to_ecliptic(t);

    // Corpo sull'eclittica a varie longitudini (lat eclittica = 0 per costruzione)
    for (double lon_deg : {0.0, 45.0, 90.0, 180.0, 270.0}) {
        const double L = lon_deg * M_PI / 180.0;
        const Vector3d in_ecl(std::cos(L), std::sin(L), 0.0);   // lat 0
        const Vector3d in_gcrf = ecl2eq * in_ecl;               // come lo da' l'effemeride

        // In GCRF la latitudine NON e' zero (c'e' l'obliquita')
        const double lat_gcrf = std::asin(in_gcrf.z() / in_gcrf.norm()) * 180.0 / M_PI;
        // La rotazione del N-body deve riportarlo a lat ~0
        const Vector3d back = eq2ecl * in_gcrf;
        const double lat_ecl = std::asin(back.z() / back.norm()) * 180.0 / M_PI;

        EXPECT_LT(std::abs(lat_ecl), 1e-9)
            << "lon=" << lon_deg << " lat_gcrf=" << lat_gcrf << " lat_ecl=" << lat_ecl;
    }
}

// Sanity sul verso: l'obliquita' recuperata deve essere ~23.44 deg, non -23 o 47.
TEST(NBodyFrame, ObliquitySignIsCorrect) {
    const auto t = time::EpochTDB::from_mjd(60000.0);
    const Matrix3d ecl2eq = coordinates::ReferenceFrame::ecliptic_to_j2000(t);
    const Vector3d pole_ecl(0, 0, 1);              // polo eclittico
    const Vector3d pole_in_eq = ecl2eq * pole_ecl;
    const double tilt = std::acos(std::clamp(pole_in_eq.z(), -1.0, 1.0)) * 180.0 / M_PI;
    EXPECT_NEAR(tilt, 23.4392911, 0.01);           // obliquita' J2000
}
