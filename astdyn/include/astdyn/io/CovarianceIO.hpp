#ifndef ASTDYN_IO_COVARIANCE_IO_HPP
#define ASTDYN_IO_COVARIANCE_IO_HPP

#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

namespace astdyn::io {

/**
 * @brief Utility for reading orbital covariance matrices.
 */
class CovarianceIO {
public:
    /**
     * @brief Read a 6x6 covariance matrix from a file.
     * Supports raw CSV (6 rows of 6 values) or OrbFit .cor format.
     */
    static Eigen::Matrix<double, 6, 6> read_file(const std::string& path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("CovarianceIO: cannot open file " + path);
        }

        Eigen::Matrix<double, 6, 6> cov = Eigen::Matrix<double, 6, 6>::Zero();
        std::string line;
        std::vector<double> values;

        while (std::getline(file, line)) {
            // Skip comments and headers
            if (line.empty() || line[0] == '!' || line[0] == '#') continue;
            if (line.find("COR") != std::string::npos) continue; // Skip AstDyS header

            std::stringstream ss(line);
            double val;
            while (ss >> val) {
                values.push_back(val);
                if (ss.peek() == ',' || ss.peek() == ';') ss.ignore();
            }
        }

        if (values.size() == 36) {
            // Full 6x6
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    cov(i, j) = values[i * 6 + j];
                }
            }
        } else if (values.size() == 21) {
            // Lower triangular (common in .cor)
            int k = 0;
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j <= i; ++j) {
                    cov(i, j) = values[k++];
                    cov(j, i) = cov(i, j);
                }
            }
        } else {
            throw std::runtime_error("CovarianceIO: expected 21 or 36 values, found " + std::to_string(values.size()));
        }

        return cov;
    }
};

} // namespace astdyn::io

#endif // ASTDYN_IO_COVARIANCE_IO_HPP
