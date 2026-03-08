#ifndef ASTDYN_ORBIT_FITTER_HPP
#define ASTDYN_ORBIT_FITTER_HPP

#include "astdyn/orbit_determination/DifferentialCorrector.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/core/physics_state.hpp"
#include <memory>

namespace astdyn::orbit_determination {

/**
 * @brief High-level service for orbit fitting.
 * 
 * Encapsulates the orchestration of ResidualCalculator, STM, and DifferentialCorrector.
 */
template <typename Frame>
class OrbitFitter {
public:
    OrbitFitter(std::shared_ptr<ephemeris::PlanetaryEphemeris> eph,
                std::shared_ptr<propagation::Propagator> prop)
        : ephemeris_(eph), propagator_(prop) 
    {
        auto residual_calc = std::make_shared<ResidualCalculator<Frame>>(ephemeris_, propagator_);
        auto stm_computer = std::make_shared<StateTransitionMatrix<Frame>>(propagator_);
        corrector_ = std::make_unique<DifferentialCorrector<Frame>>(residual_calc, stm_computer);
    }

    void set_corrections(bool aberration, bool light_time) {
        auto rc = corrector_->get_residual_calculator();
        if (rc) {
            rc->set_aberration_correction(aberration);
            rc->set_light_time_correction(light_time);
        }
    }

    DifferentialCorrectorResult<Frame> fit(
        const std::vector<observations::OpticalObservation>& obs,
        const physics::CartesianStateTyped<Frame>& initial_state,
        const DifferentialCorrectorSettings& settings) 
    {
        return corrector_->fit(obs, initial_state, settings);
    }

private:
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
    std::shared_ptr<propagation::Propagator> propagator_;
    std::unique_ptr<DifferentialCorrector<Frame>> corrector_;
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_FITTER_HPP
