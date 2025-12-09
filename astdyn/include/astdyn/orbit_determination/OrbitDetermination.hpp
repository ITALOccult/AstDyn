/**
 * @file OrbitDetermination.hpp
 * @brief Complete orbit determination system
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * Integrates all components for orbit fitting:
 * - STMPropagator
 * - ResidualCalculator  
 * - LeastSquaresFitter
 * - Parsers (EQ1, RWO)
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_HPP
#define ASTDYN_ORBIT_DETERMINATION_HPP

#include "astdyn/propagation/STMPropagator.hpp"
#include "astdyn/orbit_determination/ResidualCalculator.hpp"
#include "astdyn/orbit_determination/LeastSquaresFitter.hpp"
#include "astdyn/io/parsers/OrbFitEQ1Parser.hpp"
#include "astdyn/io/parsers/AstDysRWOParser.hpp"
#include <string>
#include <vector>

namespace astdyn::orbit_determination {

/**
 * @brief Complete orbit determination system
 * 
 * Usage:
 * ```cpp
 * OrbitDetermination od;
 * od.load_elements("asteroid.eq1");
 * od.load_observations("asteroid.rwo");
 * auto result = od.fit();
 * ```
 */
class OrbitDetermination {
public:
    /**
     * @brief Construct orbit determination system
     */
    OrbitDetermination();
    
    /**
     * @brief Load orbital elements from .eq1 file
     */
    void load_elements(const std::string& eq1_file);
    
    /**
     * @brief Load observations from .rwo file
     */
    void load_observations(const std::string& rwo_file, size_t max_obs = 0);
    
    /**
     * @brief Set initial state manually
     */
    void set_initial_state(const Eigen::Vector<double, 6>& state, double epoch_mjd);
    
    /**
     * @brief Perform orbit fit
     */
    FitResult fit();
    
    /**
     * @brief Set fitter parameters
     */
    void set_max_iterations(int max_iter) { max_iterations_ = max_iter; }
    void set_tolerance(double tol) { tolerance_ = tol; }
    void set_outlier_threshold(double sigma) { outlier_threshold_ = sigma; }
    
    /**
     * @brief Get current state
     */
    Eigen::Vector<double, 6> get_state() const { return state_; }
    double get_epoch() const { return epoch_mjd_; }
    
private:
    // State
    Eigen::Vector<double, 6> state_;
    double epoch_mjd_;
    
    // Observations (converted from OpticalObservation)
    std::vector<Observation> observations_;
    
    // Components
    std::unique_ptr<propagation::STMPropagator> stm_propagator_;
    std::unique_ptr<ResidualCalculator> residual_calculator_;
    std::unique_ptr<LeastSquaresFitter> fitter_;
    
    // Parameters
    int max_iterations_ = 10;
    double tolerance_ = 1e-6;
    double outlier_threshold_ = 3.0;
    
    /**
     * @brief Convert Keplerian elements to Cartesian state
     */
    Eigen::Vector<double, 6> elements_to_cartesian(
        const io::IOrbitParser::OrbitalElements& elements
    );
    
    /**
     * @brief Convert OpticalObservation to Observation
     */
    Observation convert_observation(
        const io::IObservationParser::OpticalObservation& opt_obs
    );
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_HPP
