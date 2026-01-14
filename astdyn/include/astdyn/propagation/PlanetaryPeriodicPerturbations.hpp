/**
 * @file PlanetaryPeriodicPerturbations.hpp
 * @brief Analytical periodic perturbations for mean to osculating conversion
 * @author ITALOccult AstDyn Team
 * @date 2025-01-14
 * 
 * This module implements the Milani-Knezevic theory (MK90/MK94) for 
 * periodic perturbations from major planets (primarily Jupiter and Saturn).
 * Ported from OrbFit src/orb9/short9.f90 and src/orb9/lcder.f90.
 */

#ifndef ASTDYN_PROPAGATION_PLANETARY_PERIODIC_PERTURBATIONS_HPP
#define ASTDYN_PROPAGATION_PLANETARY_PERIODIC_PERTURBATIONS_HPP

#include "astdyn/propagation/OrbitalElements.hpp"
#include <vector>
#include <array>

namespace astdyn::propagation {

/**
 * @brief Integer representation of a Hamiltonian term for direct perturbations
 * Corresponds to the 17-integer record in OrbFit's hamil.dir
 */
struct DirectHamTerm {
    int mul;       ///< Multiplier (ih(m,1))
    int div;       ///< Divisor (ih(m,2))
    int i_lev;     ///< Index for Leverrier coefficient (ih(m,3))
    int sign_pow;  ///< Power of -1 (ih(m,4))
    int i_pow;     ///< Power of i (ih(m,5))
    int e_pow;     ///< Power of asteroid eccentricity (ih(m,6))
    int ep_pow;    ///< Power of planetary eccentricity (ih(m,7))
    int s_pow;     ///< Power of asteroid sine inclination (ih(m,8))
    int sp_pow;    ///< Power of planetary sine inclination (ih(m,9))
    int si2_pow;   ///< Power of asteroid sin(i/2)^2 (ih(m,10))
    int si2p_pow;  ///< Power of planetary sin(i/2)^2 (ih(m,11))
    int k_p;       ///< Coefficient for planetary mean motion in div (ih(m,12))
    int k_a;       ///< Coefficient for asteroid mean motion in div (ih(m,13))
    int m1, m2, m3, m4; ///< Coefficients for longitudes in argument (ih(m,14-17))
};

/**
 * @brief Integer representation of a Hamiltonian term for indirect perturbations
 * Corresponds to the 15-integer record in OrbFit's hamil.ind
 */
struct IndirectHamTerm {
    int mul;       ///< Multiplier (ihi(m,1))
    int div;       ///< Divisor (ihi(m,2))
    int sign_pow;  ///< Power of -1 (ihi(m,3))
    int e_pow;     ///< Power of asteroid eccentricity (ihi(m,4))
    int ep_pow;    ///< Power of planetary eccentricity (ihi(m,5))
    int s_pow;     ///< Power of asteroid sine inclination (ihi(m,6))
    int sp_pow;    ///< Power of planetary sine inclination (ihi(m,7))
    int si2_pow;   ///< Power of asteroid sin(i/2)^2 (ihi(m,8))
    int si2p_pow;  ///< Power of planetary sin(i/2)^2 (ihi(m,9))
    int k_p;       ///< Coefficient for planetary mean motion (ihi(m,10))
    int k_a;       ///< Coefficient for asteroid mean motion (ihi(m,11))
    int m1, m2, m3, m4; ///< Coefficients for longitudes in argument (ihi(m,12-15))
};

/**
 * @brief Class for calculating analytical periodic perturbations
 */
class PlanetaryPeriodicPerturbations {
public:
    PlanetaryPeriodicPerturbations();

    /**
     * @brief Compute the periodic corrections for a given state
     * @param elements Mean Keplerian elements (AU, rad)
     * @param mjd Modified Julian Date (TDB)
     * @return Array of 6 corrections to Delaunay-like variables
     */
    std::array<double, 6> calculateCorrections(
        const KeplerianElements& elements,
        double mjd) const;

private:
    // Ported routines from OrbFit
    void lapco(double a, int np, double excrit, 
               std::vector<double>& b, std::vector<double>& b1, 
               std::vector<double>& b2, std::vector<double>& b3, 
               std::vector<double>& b4, std::vector<double>& b5, 
               std::vector<double>& c, std::vector<double>& c1, 
               std::vector<double>& c2, std::vector<double>& c3, 
               std::vector<double>& e, std::vector<double>& e1) const;

    void levco(int i, const std::vector<double>& b, 
               const std::vector<double>& b1, const std::vector<double>& b2, 
               const std::vector<double>& b3, const std::vector<double>& b4, 
               const std::vector<double>& c, const std::vector<double>& c1, 
               const std::vector<double>& c2, const std::vector<double>& e, 
               double a, double pai, std::vector<double>& cl) const;

    void deriv(double s, int j, double& dv, double a, double excrit, int ider) const;

    // Series data (static constants defined in CPP)
    static const std::vector<DirectHamTerm> directSeries;
    static const std::vector<IndirectHamTerm> indirectSeries;
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATION_PLANETARY_PERIODIC_PERTURBATIONS_HPP
