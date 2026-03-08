#ifndef ASTDYN_CORE_ENUMS_HPP
#define ASTDYN_CORE_ENUMS_HPP

#include <string>

namespace astdyn {

/**
 * @brief Supported Numerical Integrators
 */
enum class IntegratorType {
    RK4,    ///< Runge-Kutta 4th Order (Fixed Step)
    RKF78,  ///< Runge-Kutta-Fehlberg 7/8 (Adaptive)
    SABA4,  ///< Symplectic SABA4 (Yoshida-based)
    AAS,    ///< Adaptive Symplectic (AstDyn Precision-based)
    GAUSS,  ///< Gauss-Jackson (Legacy/Specialized)
    RADAU   ///< Radau-IIA (Implicit for stiff problems)
};

/**
 * @brief Ephemeris Data Sources
 */
enum class EphemerisType {
    Analytical, ///< Two-body or simple analytical models (Deprecated for high precision)
    DE441,      ///< JPL DE441 (Native binary provider)
    CUSTOM      ///< User-defined provider
};

/**
 * @brief Helper to convert string to IntegratorType (for IO/Config)
 */
inline IntegratorType string_to_integrator(const std::string& s) {
    if (s == "AAS") return IntegratorType::AAS;
    if (s == "RKF78") return IntegratorType::RKF78;
    if (s == "SABA4") return IntegratorType::SABA4;
    if (s == "RK4") return IntegratorType::RK4;
    if (s == "GAUSS") return IntegratorType::GAUSS;
    if (s == "RADAU") return IntegratorType::RADAU;
    return IntegratorType::RK4;
}

/**
 * @brief Helper to convert string to EphemerisType
 */
inline EphemerisType string_to_ephemeris(const std::string& s) {
    if (s == "DE441") return EphemerisType::DE441;
    if (s == "Analytical") return EphemerisType::Analytical;
    return EphemerisType::DE441;
}

} // namespace astdyn

#endif // ASTDYN_CORE_ENUMS_HPP
