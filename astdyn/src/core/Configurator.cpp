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
        if (p.contains("include_planets")) s.include_planets = p["include_planets"].get<bool>();
        if (p.contains("include_relativity")) s.include_relativity = p["include_relativity"].get<bool>();
        if (p.contains("include_moon")) s.include_moon = p["include_moon"].get<bool>();
        if (p.contains("include_asteroids")) s.include_asteroids = p["include_asteroids"].get<bool>();
        
        if (p.contains("perturb_mercury")) s.perturb_mercury = p["perturb_mercury"].get<bool>();
        if (p.contains("perturb_venus")) s.perturb_venus = p["perturb_venus"].get<bool>();
        if (p.contains("perturb_earth")) s.perturb_earth = p["perturb_earth"].get<bool>();
        if (p.contains("perturb_mars")) s.perturb_mars = p["perturb_mars"].get<bool>();
        if (p.contains("perturb_jupiter")) s.perturb_jupiter = p["perturb_jupiter"].get<bool>();
        if (p.contains("perturb_saturn")) s.perturb_saturn = p["perturb_saturn"].get<bool>();
        if (p.contains("perturb_uranus")) s.perturb_uranus = p["perturb_uranus"].get<bool>();
        if (p.contains("perturb_neptune")) s.perturb_neptune = p["perturb_neptune"].get<bool>();
        
        if (p.contains("central_body_gm")) s.central_body_gm = p["central_body_gm"].get<double>();
        
        if (p.contains("ppn_beta")) s.ppn_beta = p["ppn_beta"].get<double>();
        if (p.contains("ppn_gamma")) s.ppn_gamma = p["ppn_gamma"].get<double>();
        
        if (p.contains("asteroid_ephemeris_file")) s.asteroid_ephemeris_file = p["asteroid_ephemeris_file"].get<std::string>();
        if (p.contains("include_asteroids_list")) s.include_asteroids_list = p["include_asteroids_list"].get<std::vector<int>>();
        if (p.contains("exclude_asteroids_list")) s.exclude_asteroids_list = p["exclude_asteroids_list"].get<std::vector<int>>();
        if (p.contains("use_default_asteroid_set")) s.use_default_asteroid_set = p["use_default_asteroid_set"].get<bool>();
        if (p.contains("use_default_30_set")) s.use_default_30_set = p["use_default_30_set"].get<bool>();
        
        if (p.contains("include_earth_j2")) s.include_earth_j2 = p["include_earth_j2"].get<bool>();
        if (p.contains("include_sun_j2")) s.include_sun_j2 = p["include_sun_j2"].get<bool>();

        if (p.contains("integrate_in_ecliptic")) s.integrate_in_ecliptic = p["integrate_in_ecliptic"].get<bool>();
        if (p.contains("include_yarkovsky")) s.include_yarkovsky = p["include_yarkovsky"].get<bool>();
        if (p.contains("yarkovsky_a2")) s.yarkovsky_a2 = p["yarkovsky_a2"].get<double>();
    }
    return s;
}

orbit_determination::DifferentialCorrectorSettings Configurator::parseDifferentialCorrectorSettings(const json& j) {
    orbit_determination::DifferentialCorrectorSettings s;
    if (j.contains("differential_corrector")) {
        auto d = j["differential_corrector"];
        if (d.contains("max_iterations")) s.max_iterations = d["max_iterations"].get<int>();
        if (d.contains("convergence_tolerance_au")) s.convergence_tolerance = physics::Distance::from_au(d["convergence_tolerance_au"].get<double>());
        if (d.contains("outlier_sigma")) s.outlier_sigma = d["outlier_sigma"].get<double>();
        if (d.contains("outlier_max_sigma")) s.outlier_max_sigma = d["outlier_max_sigma"].get<double>();
        if (d.contains("outlier_min_sigma")) s.outlier_min_sigma = d["outlier_min_sigma"].get<double>();
        if (d.contains("reject_outliers")) s.reject_outliers = d["reject_outliers"].get<bool>();
        if (d.contains("compute_covariance")) s.compute_covariance = d["compute_covariance"].get<bool>();
        if (d.contains("verbose")) s.verbose = d["verbose"].get<bool>();
        if (d.contains("use_line_search")) s.use_line_search = d["use_line_search"].get<bool>();
        if (d.contains("line_search_min_alpha")) s.line_search_min_alpha = d["line_search_min_alpha"].get<double>();
        if (d.contains("rms_tolerance_arcsec")) s.rms_tolerance_arcsec = d["rms_tolerance_arcsec"].get<double>();
        if (d.contains("check_energy_barrier")) s.check_energy_barrier = d["check_energy_barrier"].get<bool>();
        if (d.contains("energy_barrier_fraction")) s.energy_barrier_fraction = d["energy_barrier_fraction"].get<double>();
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
        if (g.contains("max_iterations")) s.max_iterations = g["max_iterations"].get<int>();
        if (g.contains("tolerance_au")) s.tolerance = physics::Distance::from_au(g["tolerance_au"].get<double>());
        if (g.contains("min_separation_days")) s.min_separation = time::TimeDuration::from_days(g["min_separation_days"].get<double>());
        if (g.contains("max_separation_days")) s.max_separation = time::TimeDuration::from_days(g["max_separation_days"].get<double>());
        if (g.contains("use_light_time")) s.use_light_time = g["use_light_time"].get<bool>();
        if (g.contains("verbose")) s.verbose = g["verbose"].get<bool>();
    }
    return s;
}

close_approach::CloseApproachSettings Configurator::parseCloseApproachSettings(const json& j) {
    close_approach::CloseApproachSettings s;
    if (j.contains("close_approach")) {
        auto ca = j["close_approach"];
        if (ca.contains("detection_distance_au")) s.detection_distance = ca["detection_distance_au"].get<double>() * constants::AU;
        if (ca.contains("min_distance_au")) s.min_distance = ca["min_distance_au"].get<double>() * constants::AU;
        if (ca.contains("compute_b_plane")) s.compute_b_plane = ca["compute_b_plane"].get<bool>();
        if (ca.contains("refine_time")) s.refine_time = ca["refine_time"].get<bool>();
        if (ca.contains("time_tolerance")) s.time_tolerance = ca["time_tolerance"].get<double>();
        if (ca.contains("max_refinement_iter")) s.max_refinement_iter = ca["max_refinement_iter"].get<int>();
    }
    return s;
}

orbit_determination::GoodingIOD::Settings Configurator::parseGoodingIODSettings(const json& j) {
    orbit_determination::GoodingIOD::Settings s;
    if (j.contains("gooding_iod")) {
        auto g = j["gooding_iod"];
        if (g.contains("max_iterations")) s.max_iterations = g["max_iterations"].get<int>();
        if (g.contains("tolerance_rad")) s.tolerance_rad = g["tolerance_rad"].get<double>();
        if (g.contains("verbose")) s.verbose = g["verbose"].get<bool>();
    }
    return s;
}

orbit_determination::ExtendedKalmanFilter::Settings Configurator::parseEKFSettings(const json& j) {
    orbit_determination::ExtendedKalmanFilter::Settings s;
    if (j.contains("ekf")) {
        auto e = j["ekf"];
        if (e.contains("default_ra_sigma_arcsec")) s.default_ra_sigma = e["default_ra_sigma_arcsec"].get<double>() * constants::ARCSEC_TO_RAD;
        if (e.contains("default_dec_sigma_arcsec")) s.default_dec_sigma = e["default_dec_sigma_arcsec"].get<double>() * constants::ARCSEC_TO_RAD;
    }
    return s;
}

astrometry::AstrometricSettings Configurator::parseAstrometricSettings(const json& j) {
    astrometry::AstrometricSettings s;
    if (j.contains("astrometry")) {
        auto a = j["astrometry"];
        if (a.contains("light_time_correction")) s.light_time_correction = a["light_time_correction"];
        if (a.contains("stellar_aberration")) s.aberrazione_differenziale = a["stellar_aberration"];
        if (a.contains("frame_conversion_to_equatorial")) s.frame_conversion_to_equatorial = a["frame_conversion_to_equatorial"];
    }
    return s;
}

orbit_determination::STMSettings Configurator::parseSTMSettings(const json& j) {
    orbit_determination::STMSettings s;
    if (j.contains("stm")) {
        auto stm = j["stm"];
        if (stm.contains("use_numerical_jacobian")) s.use_numerical_jacobian = stm["use_numerical_jacobian"].get<bool>();
        if (stm.contains("differentiation_step")) s.differentiation_step = stm["differentiation_step"].get<double>();
    }
    return s;
}

void Configurator::loadFromStream(std::istream& is,
    propagation::PropagatorSettings& prop_settings,
    orbit_determination::DifferentialCorrectorSettings& dc_settings,
    orbit_determination::GaussIODSettings& iod_settings,
    orbit_determination::STMSettings& stm_settings) {
    
    json j;
    try {
        is >> j;
        prop_settings = parsePropagatorSettings(j);
        dc_settings = parseDifferentialCorrectorSettings(j);
        iod_settings = parseGaussIODSettings(j);
        stm_settings = parseSTMSettings(j);
    } catch (const json::parse_error& e) {
        std::cerr << "Configurator parse error: " << e.what() << '\n';
    }
}

void Configurator::loadFromString(const std::string& json_string,
    propagation::PropagatorSettings& prop_settings,
    orbit_determination::DifferentialCorrectorSettings& dc_settings,
    orbit_determination::GaussIODSettings& iod_settings,
    orbit_determination::STMSettings& stm_settings) {
    
    json j;
    try {
        j = json::parse(json_string);
        prop_settings = parsePropagatorSettings(j);
        dc_settings = parseDifferentialCorrectorSettings(j);
        iod_settings = parseGaussIODSettings(j);
        stm_settings = parseSTMSettings(j);
    } catch (const json::parse_error& e) {
        std::cerr << "Configurator parse error: " << e.what() << '\n';
    }
}

} // namespace astdyn::core
