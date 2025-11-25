/**
 * @file basic_usage.cpp
 * @brief Basic example showing AstDyn library usage
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 */

#include <astdyn/AstDyn.hpp>
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::constants;

void print_separator() {
    std::cout << std::string(70, '=') << "\n";
}

void demonstrate_constants() {
    print_separator();
    std::cout << "Physical and Astronomical Constants\n";
    print_separator();
    
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "\nMathematical Constants:\n";
    std::cout << "  Pi              = " << PI << "\n";
    std::cout << "  2*Pi            = " << TWO_PI << "\n";
    std::cout << "  Pi/2            = " << HALF_PI << "\n";
    std::cout << "  Deg to Rad      = " << DEG_TO_RAD << "\n";
    
    std::cout << "\nFundamental Constants:\n";
    std::cout << "  Speed of Light  = " << C_LIGHT << " km/s\n";
    std::cout << "  AU              = " << AU << " km\n";
    std::cout << "  Gauss Constant  = " << K_GAUSS << "\n";
    std::cout << "  Solar Radius    = " << R_SUN << " km\n";
    std::cout << "  Earth Radius    = " << R_EARTH << " km\n";
    
    std::cout << "\nTime Constants:\n";
    std::cout << "  JD J2000.0      = " << JD2000 << "\n";
    std::cout << "  MJD J2000.0     = " << MJD2000 << "\n";
    std::cout << "  Days per Year   = " << DAYS_PER_YEAR << "\n";
    
    std::cout << std::scientific << std::setprecision(9);
    std::cout << "\nGravitational Parameters (km³/s²):\n";
    std::cout << "  GM Sun          = " << GM_SUN << "\n";
    std::cout << "  GM Earth        = " << GM_EARTH << "\n";
    std::cout << "  GM Jupiter      = " << GM_JUPITER << "\n";
    std::cout << "  GM Mars         = " << GM_MARS << "\n";
}

void demonstrate_vectors() {
    print_separator();
    std::cout << "Vector and Matrix Operations\n";
    print_separator();
    
    // Create position vector (1 AU on x-axis)
    Vector3d position(1.0, 0.0, 0.0);
    std::cout << "\nPosition vector (AU):\n";
    std::cout << "  " << position.transpose() << "\n";
    
    // Create velocity vector
    Vector3d velocity(0.0, 0.01720209895, 0.0);  // Circular orbit velocity
    std::cout << "\nVelocity vector (AU/day):\n";
    std::cout << "  " << velocity.transpose() << "\n";
    
    // Create 6D state vector
    Vector6d state;
    state << position, velocity;
    std::cout << "\n6D State vector:\n";
    std::cout << "  Position: " << state.head<3>().transpose() << "\n";
    std::cout << "  Velocity: " << state.tail<3>().transpose() << "\n";
    
    // Vector operations
    double speed = velocity.norm();
    std::cout << "\nSpeed magnitude: " << speed << " AU/day\n";
    std::cout << "Speed in km/s:   " << speed * AU_PER_DAY_TO_KM_PER_S << " km/s\n";
    
    // Angular momentum
    Vector3d angular_momentum = position.cross(velocity);
    std::cout << "\nAngular momentum: " << angular_momentum.transpose() << "\n";
}

void demonstrate_conversions() {
    print_separator();
    std::cout << "Unit Conversions\n";
    print_separator();
    
    std::cout << std::fixed << std::setprecision(6);
    
    // Angle conversions
    double angle_deg = 45.0;
    double angle_rad = angle_deg * DEG_TO_RAD;
    std::cout << "\nAngle: " << angle_deg << " degrees = " << angle_rad << " radians\n";
    
    // Distance conversions
    double dist_au = 1.0;
    double dist_km = dist_au * AU_TO_KM;
    std::cout << "\nDistance: " << dist_au << " AU = " << dist_km << " km\n";
    
    // Velocity conversions
    double vel_au_day = 0.01720209895;
    double vel_km_s = vel_au_day * AU_PER_DAY_TO_KM_PER_S;
    std::cout << "\nVelocity: " << vel_au_day << " AU/day = " 
              << vel_km_s << " km/s\n";
    
    // Time conversions
    double jd = JD2000 + 1000.0;  // 1000 days after J2000
    double mjd = jd - 2400000.5;
    std::cout << "\nTime: JD " << jd << " = MJD " << mjd << "\n";
}

void demonstrate_result_type() {
    print_separator();
    std::cout << "Result Type for Error Handling\n";
    print_separator();
    
    // Successful result
    auto success = Result<double>::Success(42.0);
    if (success) {
        std::cout << "\nSuccess case:\n";
        std::cout << "  Value: " << success.value << "\n";
        std::cout << "  Status: " << (success.success ? "OK" : "Failed") << "\n";
    }
    
    // Failed result
    auto failure = Result<double>::Failure("Division by zero");
    if (!failure) {
        std::cout << "\nFailure case:\n";
        std::cout << "  Status: " << (failure.success ? "OK" : "Failed") << "\n";
        std::cout << "  Error: " << failure.error_message << "\n";
    }
    
    // Result with complex type
    Vector3d vec(1.0, 2.0, 3.0);
    auto vec_result = Result<Vector3d>::Success(vec);
    if (vec_result) {
        std::cout << "\nVector result:\n";
        std::cout << "  Value: " << vec_result.value.transpose() << "\n";
    }
}

void demonstrate_special_values() {
    print_separator();
    std::cout << "Special Values (NaN, Infinity)\n";
    print_separator();
    
    double normal = 1.0;
    double nan_val = NaN;
    double inf_val = Infinity;
    
    std::cout << "\nNormal value: " << normal << "\n";
    std::cout << "  isFinite: " << (isFinite(normal) ? "true" : "false") << "\n";
    std::cout << "  isNaN:    " << (isNaN(normal) ? "true" : "false") << "\n";
    
    std::cout << "\nNaN value:\n";
    std::cout << "  isFinite: " << (isFinite(nan_val) ? "true" : "false") << "\n";
    std::cout << "  isNaN:    " << (isNaN(nan_val) ? "true" : "false") << "\n";
    
    std::cout << "\nInfinity value:\n";
    std::cout << "  isFinite: " << (isFinite(inf_val) ? "true" : "false") << "\n";
    std::cout << "  isNaN:    " << (isNaN(inf_val) ? "true" : "false") << "\n";
}

int main() {
    // Initialize library
    if (!initialize()) {
        std::cerr << "ERROR: Failed to initialize AstDyn library\n";
        return 1;
    }
    
    // Print header
    std::cout << "\n";
    print_separator();
    std::cout << "AstDyn C++ - Basic Usage Example\n";
    std::cout << "Version: " << astdyn::Version::string << "\n";
    std::cout << "Build Type: " << astdyn::Config::build_type << "\n";
    std::cout << "SPICE Support: " << (astdyn::Config::use_spice ? "Enabled" : "Disabled") << "\n";
    print_separator();
    std::cout << "\n";
    
    // Run demonstrations
    demonstrate_constants();
    std::cout << "\n";
    
    demonstrate_vectors();
    std::cout << "\n";
    
    demonstrate_conversions();
    std::cout << "\n";
    
    demonstrate_result_type();
    std::cout << "\n";
    
    demonstrate_special_values();
    std::cout << "\n";
    
    // Cleanup
    shutdown();
    
    print_separator();
    std::cout << "Example completed successfully!\n";
    print_separator();
    std::cout << "\n";
    
    return 0;
}
