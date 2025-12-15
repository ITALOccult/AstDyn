// test_asteroid_fitter.cpp
// Unit tests for AsteroidFitter::fitFromConfig

#include <gtest/gtest.h>
#include "astdyn/ephemeris/AsteroidFitter.hpp"
#include "astdyn/ephemeris/AsteroidFitConfig.hpp"
#include <fstream>
#include <cstdio>

using namespace astdyn::ephemeris;

// Helper to create a temporary JSON configuration file for in‑memory test
static std::string createTempJson(const std::string& content) {
    char tmpName[] = "/tmp/asteroid_fit_config_XXXXXX.json";
    int fd = mkstemps(tmpName, 5); // keep .json suffix
    if (fd == -1) {
        throw std::runtime_error("Unable to create temporary file");
    }
    std::ofstream ofs(tmpName);
    ofs << content;
    ofs.close();
    close(fd);
    return std::string(tmpName);
}

TEST(AsteroidFitterTest, FitFromConfig_InMemorySuccess) {
    // Build a minimal JSON with orbital elements and two observation epochs
    std::string json = R"({
        "orbit": {"a":2.5,"e":0.1,"i":0.05,"Omega":1.0,"omega":0.5,"M":0.0},
        "mjd_observations": [51544.5, 51545.5],
        "outputEquatorial": true
    })";

    std::string path = createTempJson(json);
    AsteroidFitConfig cfg = loadAsteroidFitConfig(path);
    // Ensure the config was parsed correctly
    EXPECT_NEAR(cfg.orbit.a, 2.5, 1e-12);
    EXPECT_EQ(cfg.mjd_observations.size(), 2u);

    AsteroidFitResult result = AsteroidFitter::fitFromConfig(cfg);
    EXPECT_TRUE(result.success);
    // The result should contain at least the same number of positions as observations
    EXPECT_EQ(result.fitted_positions.size(), cfg.mjd_observations.size());
    // Clean up temporary file
    std::remove(path.c_str());
}

TEST(AsteroidFitterTest, FitFromConfig_NoInputFails) {
    // Empty JSON – no files, no orbit, no observations
    std::string json = "{}";
    std::string path = createTempJson(json);
    AsteroidFitConfig cfg = loadAsteroidFitConfig(path);
    AsteroidFitResult result = AsteroidFitter::fitFromConfig(cfg);
    EXPECT_FALSE(result.success);
    EXPECT_NE(result.message.find("No input data"), std::string::npos);
    std::remove(path.c_str());
}

// Note: URL‑based tests would require network access and curl; they are omitted here.
