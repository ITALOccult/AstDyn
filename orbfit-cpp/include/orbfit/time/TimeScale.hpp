/**
 * @file TimeScale.hpp
 * @brief Time scale conversions and utilities
 * @author OrbFit C++ Conversion Team
 * @date 2025-11-23
 * 
 * Conversion from Fortran time_scales.f90 module.
 * Handles conversions between different astronomical time scales:
 * UTC, UT1, TAI, TT (TDT), TDB, GPS, etc.
 */

#ifndef ORBFIT_TIME_TIMESCALE_HPP
#define ORBFIT_TIME_TIMESCALE_HPP

#include <orbfit/core/Types.hpp>
#include <orbfit/core/Constants.hpp>
#include <string>
#include <optional>

namespace orbfit {
namespace time {

// ============================================================================
// Time Scale Conversion Constants
// ============================================================================

/// TAI - UTC offset for leap seconds (changes periodically)
constexpr int LEAP_SECONDS_2025 = 37;  // As of 2025

/// TT - TAI constant offset (32.184 seconds)
constexpr double TT_MINUS_TAI = 32.184;

/// GPS - TAI constant offset (-19 seconds)
constexpr double GPS_MINUS_TAI = -19.0;

// ============================================================================
// Julian Date Conversions
// ============================================================================

/**
 * @brief Convert Modified Julian Date to Julian Date
 * 
 * @param mjd Modified Julian Date
 * @return Julian Date
 */
inline double mjd_to_jd(double mjd) {
    return mjd + 2400000.5;
}

/**
 * @brief Convert Julian Date to Modified Julian Date
 * 
 * @param jd Julian Date
 * @return Modified Julian Date
 */
inline double jd_to_mjd(double jd) {
    return jd - 2400000.5;
}

/**
 * @brief Convert MJD to calendar date (year, month, day)
 * 
 * @param mjd Modified Julian Date
 * @return Tuple of (year, month, day, fractional_day)
 */
std::tuple<int, int, int, double> mjd_to_calendar(double mjd);

/**
 * @brief Convert calendar date to MJD
 * 
 * @param year Year
 * @param month Month (1-12)
 * @param day Day of month
 * @param fraction Fractional day (0.0 - 1.0)
 * @return Modified Julian Date
 */
double calendar_to_mjd(int year, int month, int day, double fraction = 0.0);

// ============================================================================
// Time Scale Conversions
// ============================================================================

/**
 * @brief Convert UTC to TAI
 * 
 * Adds leap seconds.
 * 
 * @param mjd_utc MJD in UTC scale
 * @param leap_seconds Number of leap seconds (default: use current value)
 * @return MJD in TAI scale
 */
double utc_to_tai(double mjd_utc, int leap_seconds = LEAP_SECONDS_2025);

/**
 * @brief Convert TAI to UTC
 * 
 * Subtracts leap seconds.
 * 
 * @param mjd_tai MJD in TAI scale
 * @param leap_seconds Number of leap seconds (default: use current value)
 * @return MJD in UTC scale
 */
double tai_to_utc(double mjd_tai, int leap_seconds = LEAP_SECONDS_2025);

/**
 * @brief Convert TAI to TT (Terrestrial Time)
 * 
 * TT = TAI + 32.184 seconds
 * 
 * @param mjd_tai MJD in TAI scale
 * @return MJD in TT scale
 */
double tai_to_tt(double mjd_tai);

/**
 * @brief Convert TT to TAI
 * 
 * TAI = TT - 32.184 seconds
 * 
 * @param mjd_tt MJD in TT scale
 * @return MJD in TAI scale
 */
double tt_to_tai(double mjd_tt);

/**
 * @brief Convert UTC to TT
 * 
 * Combines UTC→TAI→TT conversions.
 * 
 * @param mjd_utc MJD in UTC scale
 * @return MJD in TT scale
 */
double utc_to_tt(double mjd_utc);

/**
 * @brief Convert TT to UTC
 * 
 * Combines TT→TAI→UTC conversions.
 * 
 * @param mjd_tt MJD in TT scale
 * @return MJD in UTC scale
 */
double tt_to_utc(double mjd_tt);

/**
 * @brief Convert TT to TDB (Barycentric Dynamical Time)
 * 
 * Uses periodic terms approximation.
 * TDB ≈ TT + 0.001658 * sin(g) + 0.000014 * sin(2*g) seconds
 * where g is the mean anomaly of Earth.
 * 
 * @param mjd_tt MJD in TT scale
 * @return MJD in TDB scale
 */
double tt_to_tdb(double mjd_tt);

/**
 * @brief Convert TDB to TT
 * 
 * Inverse of tt_to_tdb (iterative solution).
 * 
 * @param mjd_tdb MJD in TDB scale
 * @return MJD in TT scale
 */
double tdb_to_tt(double mjd_tdb);

/**
 * @brief Convert UTC to TDB
 * 
 * @param mjd_utc MJD in UTC scale
 * @return MJD in TDB scale
 */
inline double utc_to_tdb(double mjd_utc) {
    return tt_to_tdb(utc_to_tt(mjd_utc));
}

/**
 * @brief Convert TDB to UTC
 * 
 * @param mjd_tdb MJD in TDB scale
 * @return MJD in UTC scale
 */
inline double tdb_to_utc(double mjd_tdb) {
    return tt_to_utc(tdb_to_tt(mjd_tdb));
}

/**
 * @brief Convert TAI to GPS time
 * 
 * GPS = TAI - 19 seconds
 * 
 * @param mjd_tai MJD in TAI scale
 * @return MJD in GPS scale
 */
double tai_to_gps(double mjd_tai);

/**
 * @brief Convert GPS time to TAI
 * 
 * TAI = GPS + 19 seconds
 * 
 * @param mjd_gps MJD in GPS scale
 * @return MJD in TAI scale
 */
double gps_to_tai(double mjd_gps);

// ============================================================================
// UT1 Conversions (requires Earth rotation data)
// ============================================================================

/**
 * @brief Get UT1-UTC correction (DUT1)
 * 
 * Requires IERS Earth Orientation Parameters data.
 * Returns approximate value if data not loaded.
 * 
 * @param mjd_utc MJD in UTC
 * @return UT1-UTC in seconds (typically -0.9 to +0.9)
 */
double get_dut1(double mjd_utc);

/**
 * @brief Convert UTC to UT1
 * 
 * UT1 = UTC + DUT1
 * 
 * @param mjd_utc MJD in UTC scale
 * @return MJD in UT1 scale
 */
double utc_to_ut1(double mjd_utc);

/**
 * @brief Convert UT1 to UTC
 * 
 * UTC = UT1 - DUT1
 * 
 * @param mjd_ut1 MJD in UT1 scale
 * @return MJD in UTC scale
 */
double ut1_to_utc(double mjd_ut1);

// ============================================================================
// Time Utilities
// ============================================================================

/**
 * @brief Compute time difference in days
 * 
 * @param mjd1 First time (MJD)
 * @param mjd2 Second time (MJD)
 * @return Time difference in days (mjd2 - mjd1)
 */
inline double time_difference_days(double mjd1, double mjd2) {
    return mjd2 - mjd1;
}

/**
 * @brief Compute time difference in seconds
 * 
 * @param mjd1 First time (MJD)
 * @param mjd2 Second time (MJD)
 * @return Time difference in seconds (mjd2 - mjd1)
 */
inline double time_difference_seconds(double mjd1, double mjd2) {
    return (mjd2 - mjd1) * constants::SECONDS_PER_DAY;
}

/**
 * @brief Add seconds to MJD
 * 
 * @param mjd Original MJD
 * @param seconds Seconds to add
 * @return New MJD
 */
inline double add_seconds(double mjd, double seconds) {
    return mjd + seconds / constants::SECONDS_PER_DAY;
}

/**
 * @brief Add days to MJD
 * 
 * @param mjd Original MJD
 * @param days Days to add
 * @return New MJD
 */
inline double add_days(double mjd, double days) {
    return mjd + days;
}

/**
 * @brief Convert time string to MJD
 * 
 * Supported formats:
 * - "YYYY-MM-DD"
 * - "YYYY-MM-DD HH:MM:SS"
 * - "YYYY-MM-DDTHH:MM:SS" (ISO 8601)
 * 
 * @param time_string Time string
 * @param time_scale Time scale (default: UTC)
 * @return MJD, or nullopt if parsing failed
 */
std::optional<double> parse_time_string(
    const std::string& time_string,
    TimeScale time_scale = TimeScale::UTC
);

/**
 * @brief Format MJD as string
 * 
 * @param mjd Modified Julian Date
 * @param time_scale Time scale
 * @param format Format string (strftime-compatible)
 * @return Formatted time string
 */
std::string format_time(
    double mjd,
    TimeScale time_scale = TimeScale::UTC,
    const std::string& format = "%Y-%m-%d %H:%M:%S"
);

/**
 * @brief Get current system time as MJD
 * 
 * @param time_scale Desired time scale (default: UTC)
 * @return Current time as MJD
 */
double now(TimeScale time_scale = TimeScale::UTC);

/**
 * @brief Compute Julian centuries since J2000.0
 * 
 * Used for many astronomical computations.
 * T = (JD - 2451545.0) / 36525.0
 * 
 * @param jd Julian Date
 * @return Julian centuries since J2000.0
 */
inline double julian_centuries_j2000(double jd) {
    return (jd - constants::JD2000) / constants::DAYS_PER_CENTURY;
}

/**
 * @brief Compute Julian centuries since J2000.0 from MJD
 * 
 * @param mjd Modified Julian Date
 * @return Julian centuries since J2000.0
 */
inline double mjd_to_julian_centuries(double mjd) {
    return julian_centuries_j2000(mjd_to_jd(mjd));
}

// ============================================================================
// Leap Second Management
// ============================================================================

/**
 * @brief Load leap second data from file
 * 
 * @param filepath Path to leap second file
 * @return true if successful
 */
bool load_leap_seconds(const std::string& filepath);

/**
 * @brief Get number of leap seconds for a given date
 * 
 * @param mjd_utc MJD in UTC
 * @return Number of leap seconds at that date
 */
int get_leap_seconds(double mjd_utc);

/**
 * @brief Load ΔUT1 (UT1-UTC) data from IERS finals.data file
 * 
 * If filepath is empty, initializes with approximate values for 2020-2025.
 * For production use, download finals.data from:
 * https://datacenter.iers.org/data/latestVersion/finals.data
 * 
 * @param filepath Path to IERS finals.data file (empty for defaults)
 * @return true if successful
 */
bool load_dut1_data(const std::string& filepath = "");

/**
 * @brief Get ΔUT1 (UT1-UTC) for a given date
 * 
 * Uses linear interpolation between loaded values.
 * Auto-initializes with default values if not loaded.
 * 
 * @param mjd_utc MJD in UTC
 * @return ΔUT1 in seconds
 */
double get_dut1(double mjd_utc);

} // namespace time
} // namespace orbfit

#endif // ORBFIT_TIME_TIMESCALE_HPP
