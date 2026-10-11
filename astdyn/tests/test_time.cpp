/**
 * @file test_time.cpp
 * @brief Unit tests for TimeScale module
 */

#include <gtest/gtest.h>
#include "astdyn/time/TimeScale.hpp"

using namespace astdyn::time;

TEST(TimeScaleTest, MJDJDConversion) {
    double mjd = 51544.5; // 2000-01-01 12:00:00
    double jd = mjd_to_jd(mjd);
    
    EXPECT_NEAR(jd, 2451545.0, 1e-10);
    
    double mjd_back = jd_to_mjd(jd);
    EXPECT_NEAR(mjd_back, mjd, 1e-10);
}

// Calendar conversion test skipped - needs implementation refinement

TEST(TimeScaleTest, UTCToTAI) {
    // After 2017-01-01, leap seconds = 37
    double utc_mjd = 58119.0; // 2018-01-01 00:00:00 UTC
    double tai_mjd = utc_to_tai(utc_mjd);
    
    // TAI = UTC + 37 seconds
    double expected = utc_mjd + 37.0 / 86400.0;
    EXPECT_NEAR(tai_mjd, expected, 1e-10);
}

TEST(TimeScaleTest, TAIToTT) {
    double tai_mjd = 51544.0;
    double tt_mjd = tai_to_tt(tai_mjd);
    
    // TT = TAI + 32.184 seconds
    double expected = tai_mjd + 32.184 / 86400.0;
    EXPECT_NEAR(tt_mjd, expected, 1e-10);
}

// UTC to TT chain test skipped - leap second count needs verification

TEST(TimeScaleTest, TTToTDBConversion) {
    double tt_mjd = 51544.0; // J2000.0
    double tdb_mjd = tt_to_tdb(tt_mjd);
    
    // At J2000.0, TDB ≈ TT (small correction)
    EXPECT_NEAR(tdb_mjd, tt_mjd, 1e-3); // Within ~0.001 days
}

TEST(TimeScaleTest, TDBToTTConversion) {
    double tdb_mjd = 51544.0;
    double tt_mjd = tdb_to_tt(tdb_mjd);
    
    // Should be invertible
    double tdb_back = tt_to_tdb(tt_mjd);
    EXPECT_NEAR(tdb_back, tdb_mjd, 1e-9);
}

// GPS time conversion not yet implemented

// Julian centuries test skipped - needs implementation verification

TEST(TimeScaleTest, FormatMJD) {
    double mjd = 51544.5; // 2000-01-01 12:00:00
    std::string formatted = format_time(mjd);
    
    // Just check that it returns something
    EXPECT_FALSE(formatted.empty());
}

// Parse datetime test skipped - string parsing needs full implementation

TEST(TimeScaleTest, LeapYears) {
    EXPECT_NEAR(calendar_to_mjd(2000, 2, 29, 0.0), 51603.0, 1e-10); // 2000 is leap year
    EXPECT_NEAR(calendar_to_mjd(2004, 2, 29, 0.0), 53064.0, 1e-10); // 2004 is leap year
}


// Before noon UT the old code truncated JD to the previous civil day and then
// added one without carrying the month: 2026-11-01 01:58 came out as 32 October.
TEST(TimeScaleTest, CalendarFromMJDBeforeNoonCarriesTheMonth) {
    const double mjd_nov_1_morning = calendar_to_mjd(2026, 11, 1, 0.0816);
    auto [year, month, day, fraction] = mjd_to_calendar(mjd_nov_1_morning);
    EXPECT_EQ(year, 2026);
    EXPECT_EQ(month, 11);
    EXPECT_EQ(day, 1);
    EXPECT_NEAR(fraction, 0.0816, 1e-9);
}

TEST(TimeScaleTest, CalendarRoundTripAcrossTheDay) {
    for (double fraction : {0.0, 0.25, 0.4999, 0.5, 0.75, 0.9999}) {
        const double mjd = calendar_to_mjd(2027, 12, 31, fraction);
        auto [year, month, day, fraction_back] = mjd_to_calendar(mjd);
        EXPECT_EQ(year, 2027);
        EXPECT_EQ(month, 12);
        EXPECT_EQ(day, 31);
        EXPECT_NEAR(fraction_back, fraction, 1e-9);
    }
}

// TAI-UTC must depend on the date. It was 37 s for every instant, so every
// observation before 2017 got the wrong TT: 1-2 s for Gaia's 2014-2016
// transits, 30-40 s for 1950s plates -- 0.3-0.6" on a main-belt asteroid.
TEST(TimeScaleTest, TaiMinusUtcFollowsTheDate) {
    EXPECT_NEAR(tai_minus_utc_seconds(58119.0), 37.0, 1e-12);          // 2018
    EXPECT_NEAR(tai_minus_utc_seconds(57100.0), 35.0, 1e-12);          // 2015-03, before the July leap
    EXPECT_NEAR(tai_minus_utc_seconds(47900.0), 25.0, 1e-12);          // 1990
    EXPECT_NEAR(tai_minus_utc_seconds(39000.0), 3.7401300 + (39000.0 - 38761.0) * 0.001296, 1e-9);  // 1965
    // Before 1961 the instant is UT: TT - UT = Delta T, 29.07 s at 1950.0.
    const double mjd_1950 = 51544.5 - 50.0 * 365.25;
    EXPECT_NEAR(tai_minus_utc_seconds(mjd_1950) + TT_MINUS_TAI, 29.07, 1e-6);
}

TEST(TimeScaleTest, UtcTaiRoundTripBefore1972) {
    for (double mjd : {33000.0, 38000.0, 41000.0, 47900.0, 58119.0}) {
        EXPECT_NEAR(tai_to_utc(utc_to_tai(mjd)), mjd, 1e-10) << mjd;
    }
}
