/**
 * @file CovarianceIO.hpp
 * @brief Reading orbit covariance matrices, in particular AstDyS .eq1 files.
 */
#ifndef ASTDYN_IO_COVARIANCE_IO_HPP
#define ASTDYN_IO_COVARIANCE_IO_HPP

#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace astdyn::io {

/**
 * @brief An AstDyS orbit: equinoctial elements, their covariance, and the epoch.
 *
 * Units are as published by AstDyS, NOT as AstDyn's EquinoctialElements wants
 * them: a is in AU and lambda in degrees, so both the elements and the
 * covariance must be rescaled before use. See to_si().
 */
struct AstDysOrbit {
    /// a [AU], h, k, p, q [-], lambda [deg] -- ECLM J2000.
    Eigen::Matrix<double, 6, 1> elements = Eigen::Matrix<double, 6, 1>::Zero();
    /// Covariance in the same units as @c elements.
    Eigen::Matrix<double, 6, 6> covariance = Eigen::Matrix<double, 6, 6>::Zero();
    /// MJD, TDT (i.e. TT).
    double epoch_mjd_tt = 0.0;
    /// Absolute magnitude and slope, when the MAG record is present.
    double h_mag = 0.0, g_slope = 0.15;
};

/**
 * @brief Utility for reading orbital covariance matrices.
 */
class CovarianceIO {
public:
    /// AU in km, for rescaling a and its covariance.
    static constexpr double kAuKm = 149597870.700;
    static constexpr double kDegToRad = 0.017453292519943295;

    /**
     * @brief Read an AstDyS .eq1 file: elements, covariance and epoch.
     *
     * The COV records hold the 21 independent entries of the symmetric 6x6 as
     * the UPPER triangle, read row by row:
     *
     *     C00 C01 C02 C03 C04 C05 C11 C12 C13 C14 C15 C22 ... C55
     *
     * This matters: reading them as a LOWER triangle instead scrambles four of
     * the six variances and yields a matrix that is not even positive definite.
     * The (commented) RMS record in the file is the oracle -- the square roots
     * of the diagonal must reproduce it, and check_covariance() enforces that.
     */
    static AstDysOrbit read_eq1(const std::string& path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("CovarianceIO: cannot open file " + path);
        }

        AstDysOrbit orb;
        std::vector<double> cov_values;
        std::string line;
        bool have_elements = false;

        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string tag;
            if (!(ss >> tag)) continue;

            if (tag == "EQU") {
                for (int i = 0; i < 6; ++i) ss >> orb.elements(i);
                have_elements = true;
            } else if (tag == "MJD") {
                ss >> orb.epoch_mjd_tt;      // the trailing "TDT" is left in the stream
            } else if (tag == "MAG") {
                ss >> orb.h_mag >> orb.g_slope;
            } else if (tag == "COV") {
                double v;
                while (ss >> v) cov_values.push_back(v);
            }
            // RMS/EIG/WEA are commented out by AstDyS and start with '!', so the
            // tag test above skips them; NOR is the normal matrix, not needed.
        }

        if (!have_elements) {
            throw std::runtime_error("CovarianceIO: no EQU record in " + path);
        }
        // Covarianza OPZIONALE: gli .eq1 generati dal DB allnum.db (orchestratore)
        // contengono solo gli elementi, senza record COV. In quel caso la
        // covarianza resta Zero (gia' inizializzata) e il chiamante usa i soli
        // elementi (predizione geometrica, niente ellisse). Se i COV ci sono,
        // devono essere 21 (triangolo superiore della 6x6) e superare
        // check_covariance. check_covariance NON va chiamato sulla matrice Zero
        // (non definita positiva -> fallirebbe): per questo il return anticipato.
        if (cov_values.empty()) {
            return orb;   // solo elementi, covarianza Zero
        }
        if (cov_values.size() != 21) {
            throw std::runtime_error("CovarianceIO: expected 21 COV values in " + path +
                                     ", found " + std::to_string(cov_values.size()));
        }

        int k = 0;
        for (int i = 0; i < 6; ++i) {
            for (int j = i; j < 6; ++j) {          // UPPER triangle, row by row
                orb.covariance(i, j) = cov_values[k];
                orb.covariance(j, i) = cov_values[k];
                ++k;
            }
        }
        check_covariance(orb.covariance);
        return orb;
    }

    /**
     * @brief Rescale an AstDyS orbit to AstDyn's units: a in km, lambda in rad.
     *
     * The covariance transforms as C' = D C D with D = diag(AU_km, 1,1,1,1, deg2rad).
     */
    static void to_si(AstDysOrbit& orb) {
        Eigen::Matrix<double, 6, 1> d;
        d << kAuKm, 1.0, 1.0, 1.0, 1.0, kDegToRad;
        orb.elements = orb.elements.cwiseProduct(d);
        orb.covariance = d.asDiagonal() * orb.covariance * d.asDiagonal();
    }

    /**
     * @brief Reject a matrix that cannot be a covariance.
     *
     * A covariance is symmetric and positive definite by construction. Checking
     * it costs nothing and catches exactly the failure mode that a wrong
     * triangle ordering produces, which is otherwise completely silent: both
     * matrices are 6x6 and every subsequent multiplication succeeds.
     */
    static void check_covariance(const Eigen::Matrix<double, 6, 6>& C) {
        if (!C.isApprox(C.transpose(), 1e-12)) {
            throw std::runtime_error("CovarianceIO: matrix is not symmetric");
        }
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> es(C);
        if (es.eigenvalues().minCoeff() <= 0.0) {
            throw std::runtime_error(
                "CovarianceIO: matrix is not positive definite (smallest eigenvalue "
                + std::to_string(es.eigenvalues().minCoeff()) +
                "); the triangle ordering of the input is the usual cause");
        }
    }

    /**
     * @brief Read a bare 6x6 covariance: CSV (36 values) or 21 upper-triangle values.
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
            if (line.empty() || line[0] == '!' || line[0] == '#') continue;
            std::stringstream ss(line);
            std::string first;
            std::streampos save = ss.tellg();
            if (ss >> first && (first == "COV" || first == "COR")) {
                // consume the tag
            } else {
                ss.clear();
                ss.seekg(save);
            }
            double val;
            while (ss >> val) {
                values.push_back(val);
                if (ss.peek() == ',' || ss.peek() == ';') ss.ignore();
            }
        }

        if (values.size() == 36) {
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) cov(i, j) = values[i * 6 + j];
        } else if (values.size() == 21) {
            int k = 0;
            for (int i = 0; i < 6; ++i) {
                for (int j = i; j < 6; ++j) {      // UPPER triangle, as AstDyS writes it
                    cov(i, j) = values[k];
                    cov(j, i) = values[k];
                    ++k;
                }
            }
        } else {
            throw std::runtime_error("CovarianceIO: expected 21 or 36 values in " + path +
                                     ", found " + std::to_string(values.size()));
        }
        check_covariance(cov);
        return cov;
    }
};

} // namespace astdyn::io

#endif // ASTDYN_IO_COVARIANCE_IO_HPP
