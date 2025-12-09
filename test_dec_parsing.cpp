/**
 * @file test_dec_parsing.cpp
 * @brief Quick test for Dec parsing
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "Testing Dec parsing from 17030_astdys.rwo\n\n";
    
    std::ifstream file("17030_astdys.rwo");
    std::string line;
    int count = 0;
    
    while (std::getline(file, line) && count < 10) {
        if (line.empty() || line[0] == '!' || line.find("END_OF_HEADER") != std::string::npos) {
            continue;
        }
        
        std::istringstream iss(line);
        std::string fields[25];
        for (int i = 0; i < 25 && iss >> fields[i]; ++i);
        
        // RA: fields 7, 8, 9
        // Dec: fields 15, 16, 17
        
        if (fields[15].length() > 0) {
            try {
                // RA
                int ra_h = std::stoi(fields[7]);
                int ra_m = std::stoi(fields[8]);
                double ra_s = std::stod(fields[9]);
                double ra_deg = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0;
                
                // Dec
                std::string dec_d_str = fields[15];
                int dec_m = std::stoi(fields[16]);
                double dec_s = std::stod(fields[17]);
                
                char sign = dec_d_str[0];
                int dec_d = std::abs(std::stoi(dec_d_str));
                
                double dec_deg = dec_d + dec_m/60.0 + dec_s/3600.0;
                if (sign == '-') dec_deg = -dec_deg;
                
                std::cout << std::fixed << std::setprecision(6);
                std::cout << "RA: " << std::setw(10) << ra_deg 
                          << "  Dec: " << std::setw(10) << dec_deg << "\n";
                
                count++;
            } catch (...) {
                continue;
            }
        }
    }
    
    return 0;
}
