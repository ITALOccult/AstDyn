/**
 * @file test_eq1_parser.cpp
 * @brief Test del parser EQ1
 */

#include <astdyn/io/EQ1Parser.hpp>
#include <iostream>
#include <iomanip>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <file.eq1>\n";
        return 1;
    }

    try {
        std::cout << "Parsing " << argv[1] << "...\n\n";

        auto elem = astdyn::io::EQ1Parser::parse(argv[1]);

        std::cout << "Object: " << elem.object_name << "\n";
        std::cout << "Epoch: MJD " << std::fixed << std::setprecision(6) << elem.epoch_mjd_tdb << " TDB\n\n";

        std::cout << "Equinoctial Elements:\n";
        std::cout << "  a      = " << std::setprecision(10) << elem.a << " AU\n";
        std::cout << "  h      = " << elem.h << "\n";
        std::cout << "  k      = " << elem.k << "\n";
        std::cout << "  p      = " << elem.p << "\n";
        std::cout << "  q      = " << elem.q << "\n";
        std::cout << "  lambda = " << elem.lambda << "°\n\n";

        // Convert to Keplerian
        double a, e, i, Omega, omega, M;
        astdyn::io::EQ1Parser::equinoctial_to_keplerian(elem, a, e, i, Omega, omega, M);

        const double RAD_TO_DEG = 180.0 / 3.14159265358979323846;

        std::cout << "Keplerian Elements:\n";
        std::cout << "  a = " << std::setprecision(10) << a << " AU\n";
        std::cout << "  e = " << std::setprecision(8) << e << "\n";
        std::cout << "  i = " << std::setprecision(6) << i * RAD_TO_DEG << "°\n";
        std::cout << "  Ω = " << Omega * RAD_TO_DEG << "°\n";
        std::cout << "  ω = " << omega * RAD_TO_DEG << "°\n";
        std::cout << "  M = " << M * RAD_TO_DEG << "°\n\n";

        if (elem.magnitude != 0.0) {
            std::cout << "Magnitude: H = " << elem.magnitude << ", G = " << elem.mag_slope << "\n";
        }

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
