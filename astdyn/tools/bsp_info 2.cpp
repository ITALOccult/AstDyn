#include <astdyn/io/SPKReader.hpp>
#include <iostream>
#include <iomanip>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2) { std::cout << "Usage: bsp_info <file.bsp>\n"; return 1; }
    try {
        astdyn::io::SPKReader spk(argv[1]);
        auto segments = spk.getSegments();
        std::cout << "\n=== BSP Information Report ===\n";
        std::cout << "File: " << argv[1] << "\n";
        std::cout << "Number of segments: " << segments.size() << "\n\n";
        
        std::cout << std::left << std::setw(12) << "Body ID" 
                  << std::setw(12) << "Center ID" 
                  << std::setw(10) << "Frame" 
                  << std::setw(8) << "Type" 
                  << "Range (ET s from J2000)\n";
        std::cout << std::string(80, '-') << "\n";
        
        for (const auto& s : segments) {
            std::cout << std::left << std::setw(12) << s.body_id 
                      << std::setw(12) << s.center_id 
                      << std::setw(10) << s.frame_id 
                      << std::setw(8) << s.type 
                      << std::fixed << std::setprecision(1) << s.start_et << " to " << s.end_et << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
    return 0;
}
