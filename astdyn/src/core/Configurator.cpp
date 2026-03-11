#include "astdyn/core/Configurator.hpp"
#include "astdyn/propagation/OrbFitIntegrator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/AASIntegrator.hpp"
#include "astdyn/propagation/RadauIntegrator.hpp"
#include "astdyn/propagation/saba4_integrator.hpp"
#include "astdyn/propagation/GaussIntegrator.hpp"
#include <iostream>

namespace astdyn::core {

using json = nlohmann::json;

propagation::PropagatorSettings Configurator::parsePropagatorSettings(const json& j) {
    propagation::PropagatorSettings s;
    if (j.contains("propagator")) {
        auto p = j["propagator"];
        if (p.contains("include_planets")) s.include_planets = p["include_planets"];
        if (p.contains("include_relativity")) s.include_relativity = p["include_relativity"];
        if (p.contains("include_moon")) s.include_moon = p["include_moon"];
        if (p.contains("include_asteroids")) s.include_asteroids = p["include_asteroids"];
        
        if (p.contains("perturb_mercury")) s.perturb_mercury = p["perturb_mercury"];
        if (p.contains("perturb_venus")) s.perturb_venus = p["perturb_venus"];
        if (p.contains("perturb_earth")) s.perturb_earth = p["perturb_earth"];
        if (p.contains("perturb_mars")) s.perturb_mars = p["perturb_mars"];
        if (p.contains("perturb_jupiter")) s.perturb_jupiter = p["perturb_jupiter"];
        if (p.contains("perturb_saturn")) s.perturb_saturn = p["perturb_saturn"];
        if (p.contains("perturb_uranus")) s.perturb_uranus = p["perturb_uranus"];
        if (p.contains("perturb_neptune")) s.perturb_neptune = p["perturb_neptune"];
        
        if (p.contains("central_body_gm")) s.central_body_gm = p["central_body_gm"];
        
        if (p.contains("ppn_beta")) s.ppn_beta = p["ppn_beta"];
        if (p.contains("ppn_gamma")) s.ppn_gamma = p["ppn_gamma"];
        
        if (p.contains("asteroid_ephemeris_file")) s.asteroid_ephemeris_file = p["asteroid_ephemeris_file"];
        if (p.contains("integrate_in_ecliptic")) s.integrate_in_ecliptic = p["integrate_in_ecliptic"];
        if (p.contains("include_yarkovsky")) s.include_yarkovsky = p["include_yarkovsky"];
        if (p.contains("yarkovsky_a2")) s.yarkovsky_a2 = p["yarkovsky_a2"];
    }
    return s;
}

orbit_determination::DifferentialCorrectorSettings Configurator::parseDifferentialCorrectorSettings(const json& j) {
    orbit_determination::DifferentialCorrectorSettings s;
    if (j.contains("differential_corrector")) {
        auto d = j["differential_corrector"];
        if (d.contains("max_iterations")) s.max_iterations = d["max_iterations"];
        if (d.contains("convergence_tolerance_au")) s.convergence_tolerance = physics::Distance::from_au(d["convergence_tolerance_au"]);
        if (d.contains("outlier_sigma")) s.outlier_sigma = d["outlier_sigma"];
        if (d.contains("outlier_max_sigma")) s.outlier_max_sigma = d["outlier_max_sigma"];
        if (d.contains("outlier_min_sigma")) s.outlier_min_sigma = d["outlier_min_sigma"];
        if (d.contains("reject_outliers")) s.reject_outliers = d["reject_outliers"];
        if (d.contains("compute_covariance")) s.compute_covariance = d["compute_covariance"];
        if (d.contains("verbose")) s.verbose = d["verbose"];
        if (d.contains("use_line_search")) s.use_line_search = d["use_line_search"];
        if (d.contains("line_search_min_alpha")) s.line_search_min_alpha = d["line_search_min_alpha"];
        if (d.contains("rms_tolerance_arcsec")) s.rms_tolerance_arcsec = d["rms_tolerance_arcsec"];
        if (d.contains("check_energy_barrier")) s.check_energy_barrier = d["check_energy_barrier"];
        if (d.contains("energy_barrier_fraction")) s.energy_barrier_fraction = d["energy_barrier_fraction"];
    }
    return s;
}

std::shared_ptr<propagation::Integrator> Configurator::createIntegrator(const json& j) {
    if (j.contains("integrator")) {
        auto type = j["integrator"].value("type", "RKF78");
        double initial_step = j["integrator"].value("initial_step", 0.1);
        double tolerance = j["integrator"].value("tolerance", 1e-12);
        
        if (type == "OrbFitDP") {
            double min_step = j["integrator"].value("min_step", 1e-6);
            double max_step = j["integrator"].value("max_step", 5.0);
            return std::make_shared<propagation::OrbFitDPIntegrator>(initial_step, tolerance, min_step, max_step);
        } else if (type == "RK4") {
            return std::make_shared<propagation::RK4Integrator>(initial_step);
        } else if (type == "OrbFitRK4") {
            return std::make_shared<propagation::OrbFitRK4Integrator>(initial_step);
        } else if (type == "AAS") {
            double precision = j["integrator"].value("precision", 1e-4);
            return std::make_shared<propagation::AASIntegrator>(precision);
        } else if (type == "Radau" || type == "Radau15") {
            double min_step = j["integrator"].value("min_step", 1e-8);
            double max_step = j["integrator"].value("max_step", 100.0);
            int max_newton_iter = j["integrator"].value("max_newton_iter", 7);
            return std::make_shared<propagation::RadauIntegrator>(initial_step, tolerance, min_step, max_step, max_newton_iter);
        } else if (type == "SABA4") {
            double min_step = j["integrator"].value("min_step", 1e-6);
            double max_step = j["integrator"].value("max_step", 100.0);
            return std::make_shared<propagation::SABA4Integrator>(initial_step, tolerance, min_step, max_step);
        } else if (type == "Gauss") {
            double req_accuracy = j["integrator"].value("req_accuracy", 1e-12);
            return std::make_shared<propagation::GaussIntegrator>(req_accuracy, initial_step);
        } else {
            // Default to RKF78
            double min_step = j["integrator"].value("min_step", 1e-6);
            double max_step = j["integrator"].value("max_step", 100.0);
            return std::make_shared<propagation::RKF78Integrator>(initial_step, tolerance, min_step, max_step);
        }
    }
    // Fallback default
    return std::make_shared<propagation::RKF78Integrator>(0.1, 1e-12);
}

std::shared_ptr<propagation::Propagator> Configurator::createPropagator(const json& j, std::shared_ptr<ephemeris::PlanetaryEphemeris> ephem) {
    auto integr = createIntegrator(j);
    auto settings = parsePropagatorSettings(j);
    return std::make_shared<propagation::Propagator>(integr, ephem, settings);
}

orbit_determination::GaussIODSettings Configurator::parseGaussIODSettings(const json& j) {
    orbit_determination::GaussIODSettings s;
    if (j.contains("gauss_iod")) {
        auto g = j["gauss_iod"];
        if (g.contains("max_iterations")) s.max_iterations = g["max_iterations"];
        if (g.contains("tolerance_au")) s.tolerance = physics::Distance::from_au(g["tolerance_au"]);
        if (g.contains("min_separation_days")) s.min_separation_days = g["min_separation_days"];
        if (g.contains("max_separation_days")) s.max_separation_days = g["max_separation_days"];
        if (g.contains("use_light_time")) s.use_light_time = g["use_light_time"];
        if (g.contains("verbose")) s.verbose = g["verbose"];
    }
    return s;
}

void Configurator::loadFromStream(std::istream& is,
    propagation::PropagatorSettings& prop_settings,
    orbit_determination::DifferentialCorrectorSettings& dc_settings,
    orbit_determination::GaussIODSettings& iod_settings) {
    
    json j;
    try {
        is >> j;
        prop_settings = parsePropagatorSettings(j);
        dc_settings = parseDifferentialCorrectorSettings(j);
        iod_settings = parseGaussIODSettings(j);
    } catch (const json::parse_error& e) {
        std::cerr << "Configurator parse error: " << e.what() << '\n';
    }
}

void Configurator::loadFromString(const std::string& json_string,
    propagation::PropagatorSettings& prop_settings,
    orbit_determination::DifferentialCorrectorSettings& dc_settings,
    orbit_determination::GaussIODSettings& iod_settings) {
    
    json j;
    try {
        j = json::parse(json_string);
        prop_settings = parsePropagatorSettings(j);
        dc_settings = parseDifferentialCorrectorSettings(j);
        iod_settings = parseGaussIODSettings(j);
    } catch (const json::parse_error& e) {
        std::cerr << "Configurator parse error: " << e.what() << '\n';
    }
}

} // namespace astdyn::core
