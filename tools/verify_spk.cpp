#include "astdyn/io/SPKReader.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <filesystem>

// Simple diagnostic tool to inspecting SPK files using native reader
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./verify_spk <bsp_file_path>\n";
        return 1;
    }
    
    std::string path = argv[1];
    std::cout << "Inspecting SPK file: " << path << "\n";
    
    try {
        astdyn::io::SPKReader reader(path);
        
        // Access protected/private members? No, SPKReader encapsulates them.
        // But we added debug prints in the previous step.
        // Let's try to query some standard bodies to see if they exist.
        
        std::vector<int> bodies_to_check = {
            1, // Mercury Bary
            2, // Venus Bary
            3, // Earth-Moon Bary
            4, // Mars Bary
            399, // Earth
            499, // Mars
            10, // Sun
        };
        
        for (int id : bodies_to_check) {
            std::cout << "Checking Body " << id << "... ";
            try {
                // Try to get state at J2000 (0 ET)
                // If the file is Part 2 (1969+), this should work.
                // If Part 1 (-13000 to 1969), J2000 is OUT of range.
                // So let's try an epoch that should be in Part 2: 2026 (approx 8.2e8 sec)
                
                double et_test = 5.59e8; // Approx 2017 (start of integration)
                auto state = reader.getState(id, et_test);
                std::cout << "FOUND. State(5.59e8): " << state.head<3>().transpose() << "\n";
            } catch (const std::exception& e) {
                std::cout << "FAILED: " << e.what() << "\n";
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "CRITICAL ERROR loading file: " << e.what() << "\n";
    }
    
    return 0;
}
