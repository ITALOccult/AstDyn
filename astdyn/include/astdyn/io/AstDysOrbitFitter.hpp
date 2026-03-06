#ifndef ASTDYN_IO_ASTDYSORBITFITTER_HPP
#define ASTDYN_IO_ASTDYSORBITFITTER_HPP

#include <string>
#include <vector>
#include <optional>
#include "astdyn/core/physics_state.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/AstDynEngine.hpp"

namespace astdyn {
namespace io {

struct AstDysFitResult {
    physics::KeplerianStateTyped<core::ECLIPJ2000> fitted_orbit;
    bool converged = false;
    int num_iterations = 0;
    double rms_ra = 0.0;
    double rms_dec = 0.0;
    double chi_squared = 0.0;
    int num_observations_loaded = 0;
    int num_observations_used = 0;
    int num_outliers = 0;
    
    std::optional<physics::KeplerianStateTyped<core::ECLIPJ2000>> reference_orbit;
    std::optional<double> delta_a_km;
    std::optional<double> delta_e;
    std::optional<double> delta_i_arcsec;
};

class AstDysOrbitFitter {
public:
    AstDysOrbitFitter();
    
    void set_observations_file(const std::string& filename, const std::string& format = "auto");
    void set_observations(const std::vector<observations::OpticalObservation>& observations);
    
    void set_elements_file(const std::string& filename, const std::string& format = "auto");
    void set_elements(const physics::KeplerianStateTyped<core::ECLIPJ2000>& elements);
    
    void set_config_file(const std::string& filename);
    void set_reference_orbit(const physics::KeplerianStateTyped<core::ECLIPJ2000>& elements);
    
    void set_verbose(bool verbose) { verbose_ = verbose; }
    
    AstDysFitResult fit();
    
    static AstDysFitResult fit_from_astdys(const std::string& object_name, 
                                          const std::string& download_dir = ".");

private:
    std::vector<observations::OpticalObservation> observations_;
    std::optional<physics::KeplerianStateTyped<core::ECLIPJ2000>> initial_elements_;
    std::optional<physics::KeplerianStateTyped<core::ECLIPJ2000>> reference_orbit_;
    std::optional<AstDynConfig> config_;
    bool verbose_;
    
    void load_rwo_file(const std::string& filename);
    void load_mpc_file(const std::string& filename);
    void load_eq1_file(const std::string& filename);
    void load_oel_file(const std::string& filename);
    void load_json_file(const std::string& filename);
    
    std::string detect_format(const std::string& filename);
    
    static physics::KeplerianStateTyped<core::ECLIPJ2000> equinoctial_to_keplerian(
        double a, double h, double k, double p, double q, double lambda, double mjd);
};

} // namespace io
} // namespace astdyn

#endif // ASTDYN_IO_ASTDYSORBITFITTER_HPP
