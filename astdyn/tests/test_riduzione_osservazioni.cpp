/**
 * @file test_riduzione_osservazioni.cpp
 * @brief Inputs of the orbit fit: observatory sites and the AstDyS .rwo reader.
 *
 * Each case is a defect measured on (3830) against the residuals AstDyS
 * publishes in the same .rwo file:
 *  - sites built from the MPC parallax constants, without a second flattening
 *    (~15 mas, sign of the hemisphere);
 *  - the FCCT14 catalogue debiasing subtracted (0.10" -> 0.05" of AstDyS);
 *  - satellite observations (TESS, WISE) read with their observer position.
 */

#include <gtest/gtest.h>
#include "astdyn/observations/ObservatoryDatabase.hpp"
#include "astdyn/observations/RWOReader.hpp"
#include "astdyn/observations/NightlyDeweighting.hpp"
#include "astdyn/observations/GaiaSSOReader.hpp"
#include "astdyn/orbit_determination/StationWeights.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

using namespace astdyn;

namespace {

constexpr double WGS84_A_M = 6378137.0;

const char* kLineaConBias =
    " 3830      O C   1998 12 16.39428       1.000E-05 07 52 24.900  1.500E-01    1.022 F   -0.264"
    "    1.346 +14 18 39.20  1.000E-01    0.883 F   -0.670   -0.123 16.6 R  0.70   0.02   a 699      1.33 1 1";

const char* kSatellite =
    " 3830      S S   2010 05 01.11811       1.000E-05 09 00 17.360  1.500E-01    1.000 F    0.118"
    "   -0.049 +06 30 15.20  1.000E-01    1.000 F    0.090   -0.051                       L C51      0.07 1 0\n"
    " 3830      S s   2010 05 01.11811 1      -3883.985400000000       5589.788600000001       1169.687200000000 C51\n";

double hms_to_rad(int h, int m, double s) { return (h + m / 60.0 + s / 3600.0) * 15.0 * constants::DEG_TO_RAD; }
double dms_to_rad(int d, int m, double s) { return (d + m / 60.0 + s / 3600.0) * constants::DEG_TO_RAD; }

}  // namespace

TEST(RiduzioneOsservazioni, SitoDalleCostantiDiParallasseSenzaDoppioSchiacciamento) {
    const std::string path = ::testing::TempDir() + "obscodes_g96.txt";
    std::ofstream(path) << "G96 249.211280.845107+0.533611University of Arizona Mt. Lemmon Survey\n";
    auto& db = observations::ObservatoryDatabase::getInstance();
    db.loadFromMPCFile(path);
    const auto site = db.getObservatory("G96");
    ASSERT_TRUE(site.has_value());

    const double lon_rad = 249.21128 * constants::DEG_TO_RAD;
    EXPECT_NEAR(site->position.x_si(), WGS84_A_M * 0.845107 * std::cos(lon_rad), 1.0);
    EXPECT_NEAR(site->position.y_si(), WGS84_A_M * 0.845107 * std::sin(lon_rad), 1.0);
    EXPECT_NEAR(site->position.z_si(), WGS84_A_M * 0.533611, 1.0);
    std::remove(path.c_str());
}

TEST(RiduzioneOsservazioni, DebiasingDelCatalogoSottratto) {
    std::istringstream in(std::string(kLineaConBias) + "\n");
    const auto obs = observations::RWOReader::readStream(in);
    ASSERT_EQ(obs.size(), 1u);

    // The Bias columns are in arcsec of great circle and are SUBTRACTED.
    const double dec_raw = dms_to_rad(14, 18, 39.20);
    const double expected_dec = dec_raw + 0.670 * constants::ARCSEC_TO_RAD;
    const double expected_ra = hms_to_rad(7, 52, 24.900) + 0.264 * constants::ARCSEC_TO_RAD / std::cos(expected_dec);
    EXPECT_NEAR(obs[0].dec.to_rad(), expected_dec, 1e-4 * constants::ARCSEC_TO_RAD);
    EXPECT_NEAR(obs[0].ra.to_rad(), expected_ra, 1e-3 * constants::ARCSEC_TO_RAD);
    EXPECT_EQ(obs[0].observatory_code, "699");
}

TEST(RiduzioneOsservazioni, OsservazioneDaSatelliteConPosizioneDellOsservatore) {
    std::istringstream in(kSatellite);
    const auto obs = observations::RWOReader::readStream(in);
    ASSERT_EQ(obs.size(), 1u);
    EXPECT_EQ(obs[0].observatory_code, "C51");
    ASSERT_TRUE(obs[0].space_observer_gcrf.has_value());
    EXPECT_NEAR(obs[0].space_observer_gcrf->x_si(), -3883985.4, 1e-3);
    EXPECT_NEAR(obs[0].space_observer_gcrf->y_si(), 5589788.6, 1e-3);
    EXPECT_NEAR(obs[0].space_observer_gcrf->z_si(), 1169687.2, 1e-3);
}

TEST(RiduzioneOsservazioni, RigaSatelliteSenzaPosizioneScartata) {
    std::istringstream in(std::string(kSatellite).substr(0, std::string(kSatellite).find('\n') + 1));
    EXPECT_TRUE(observations::RWOReader::readStream(in).empty());
}

TEST(RiduzioneOsservazioni, PesiPerNotteVeres2017) {
    std::vector<observations::OpticalObservation> obs;
    auto make = [](const std::string& code, double mjd) {
        observations::OpticalObservation o;
        o.observatory_code = code;
        o.time = time::EpochUTC::from_mjd(mjd);
        o.sigma_ra = astrometry::Angle::from_arcsec(0.6);
        o.sigma_dec = astrometry::Angle::from_arcsec(0.6);
        return o;
    };
    for (int k = 0; k < 16; ++k) obs.push_back(make("C57", 59282.1 + k * 0.01));   // one night, 16 obs
    for (int k = 0; k < 4; ++k) obs.push_back(make("T05", 59282.3 + k * 0.01));    // 4: untouched
    obs.push_back(make("C57", 59283.2));                                            // next night, alone

    EXPECT_EQ(observations::deweight_nightly_batches(obs), 16);
    EXPECT_NEAR(obs[0].sigma_ra.to_arcsec(), 0.6 * 2.0, 1e-12);   // sqrt(16/4) = 2
    EXPECT_NEAR(obs[16].sigma_ra.to_arcsec(), 0.6, 1e-12);
    EXPECT_NEAR(obs[20].sigma_dec.to_arcsec(), 0.6, 1e-12);
}

TEST(RiduzioneOsservazioni, GaiaFprCovarianzaCorrelataEInflataPerTransito) {
    const std::string path = ::testing::TempDir() + "gaia_fpr.csv";
    std::ofstream(path)
        << "transit_id,epoch_utc,ra,dec,ra_error_systematic,dec_error_systematic,ra_dec_correlation_systematic,"
           "ra_error_random,dec_error_random,ra_dec_correlation_random,x_gaia_geocentric,y_gaia_geocentric,"
           "z_gaia_geocentric,is_rejected\n"
        << "7,1721.5,119.2,17.0,2.0,1.0,0.0,400.0,300.0,-0.99,0.0092,-0.0023,0.0002,false\n"
        << "7,1721.6,119.2,17.0,2.0,1.0,0.0,400.0,300.0,-0.99,0.0092,-0.0023,0.0002,false\n"
        << "8,1722.0,119.3,17.1,2.0,1.0,0.0,400.0,300.0,-0.99,0.0092,-0.0023,0.0002,true\n";
    const auto obs = observations::read_gaia_sso_csv(path);
    ASSERT_EQ(obs.size(), 2u);                                   // rejected record skipped
    ASSERT_TRUE(obs[0].radec_covariance_arcsec2.has_value());
    const auto& c = *obs[0].radec_covariance_arcsec2;
    // systematic counted once per transit: 2 observations -> variance x 2
    EXPECT_NEAR(c(0, 0), 2.0 * 0.002 * 0.002 + 0.4 * 0.4, 1e-12);
    EXPECT_NEAR(c(0, 1), -0.99 * 0.4 * 0.3, 1e-12);
    EXPECT_EQ(obs[0].observatory_code, "258");
    EXPECT_NEAR(obs[0].time.mjd(), 2455197.5 + 1721.5 - 2400000.5, 1e-9);
    ASSERT_TRUE(obs[0].space_observer_gcrf.has_value());
    EXPECT_NEAR(obs[0].space_observer_gcrf->x_si(), 0.0092 * constants::AU * 1000.0, 1.0);
    std::remove(path.c_str());
}

TEST(RiduzioneOsservazioni, SatelliteConNotaInColonna15) {
    // TESS 2023: the note "0" in column 15 must not shift the fields.
    std::istringstream in(
        " 316       S S 0 2023 12 05.005868      1.000E-06 11 32 16.824  1.500E-02    0.600 F    0.000    0.155 "
        "+04 14 40.88  1.000E-02    0.600 F    0.000    0.386 15.6 G  0.70  -0.09   V C57      0.69 1 1\n"
        " 316       S s 0 2023 12 05.00586 1     -22382.197000000004     220051.081000000006      97190.272999999986 C57\n");
    const auto obs = observations::RWOReader::readStream(in);
    ASSERT_EQ(obs.size(), 1u);
    ASSERT_TRUE(obs[0].space_observer_gcrf.has_value());
    EXPECT_NEAR(obs[0].space_observer_gcrf->x_si(), -22382197.0, 1.0);
    EXPECT_NEAR(obs[0].space_observer_gcrf->y_si(), 220051081.0, 1.0);
    EXPECT_NEAR(obs[0].space_observer_gcrf->z_si(), 97190273.0, 1.0);
}

TEST(RiduzioneOsservazioni, PesiPerStazioneDallaDispersioneDeiResidui) {
    std::vector<orbit_determination::ObservationResidual> res;
    auto add = [&res](const std::string& code, double z) {
        orbit_determination::ObservationResidual r;
        r.observatory_code = code; r.outlier = false;
        r.normalized_ra = z; r.normalized_dec = -z;
        res.push_back(r);
    };
    // "A": residuals twice their sigma -> factor 2; "B": too few -> untouched
    for (int k = 0; k < 200; ++k) add("A", 2.0 * (k % 2 ? 0.6745 : -0.6745));
    for (int k = 0; k < 5; ++k) add("B", 9.0);
    const auto factors = orbit_determination::station_sigma_factors(res);
    ASSERT_EQ(factors.count("A"), 1u);
    EXPECT_EQ(factors.count("B"), 0u);
    EXPECT_NEAR(factors.at("A"), 2.0, 0.01);

    std::vector<observations::OpticalObservation> obs(2);
    obs[0].observatory_code = "A"; obs[0].sigma_ra = astrometry::Angle::from_arcsec(0.5);
    obs[1].observatory_code = "A"; obs[1].radec_covariance_arcsec2 = Eigen::Matrix2d::Identity();
    obs[1].sigma_ra = astrometry::Angle::from_arcsec(0.5);
    EXPECT_EQ(orbit_determination::apply_station_factors(obs, factors), 1);   // Gaia-like untouched
    EXPECT_NEAR(obs[0].sigma_ra.to_arcsec(), 1.0, 0.01);
    EXPECT_NEAR(obs[1].sigma_ra.to_arcsec(), 0.5, 1e-12);
}
