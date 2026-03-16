#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "astdyn/AstDyn.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/io/MPCParser.hpp"
#include "astdyn/orbit_determination/GoodingIOD.hpp"
#include "astdyn/orbit_determination/ExtendedKalmanFilter.hpp"
#include "astdyn/orbit_determination/ResidualAnalysis.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include "astdyn/orbit_determination/CovariancePropagator.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/CelestialBody.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/orbit_determination/DifferentialCorrector.hpp"
#include "astdyn/api/OrbitFitAPI.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/astrometry/Astrometry.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"

namespace py = pybind11;
using namespace astdyn;

PYBIND11_MODULE(pyastdyn, m) {
    m.doc() = "Python bindings for AstDyn 3.0 - High-precision Asteroid Dynamics";

    // --- Time Module ---
    py::class_<time::EpochTDB>(m, "EpochTDB")
        .def_static("from_jd", &time::EpochTDB::from_jd)
        .def_static("from_mjd", &time::EpochTDB::from_mjd)
        .def_static("now", []() { return time::to_tdb(time::EpochTT::from_mjd(time::now(TimeScale::TT))); })
        .def("jd", &time::EpochTDB::jd)
        .def("mjd", &time::EpochTDB::mjd)
        .def("__repr__", [](const time::EpochTDB& self) {
            return "<pyastdyn.EpochTDB MJD=" + std::to_string(self.mjd()) + ">";
        });

    py::class_<time::EpochUTC>(m, "EpochUTC")
        .def_static("from_mjd", &time::EpochUTC::from_mjd)
        .def("mjd", &time::EpochUTC::mjd)
        .def("to_tdb", [](const time::EpochUTC& self) { return time::to_tdb(self); });

    // --- Astrometry Module ---
    py::class_<astrometry::Angle>(m, "Angle")
        .def_static("from_deg", &astrometry::Angle::from_deg)
        .def_static("from_rad", &astrometry::Angle::from_rad)
        .def_static("from_arcsec", &astrometry::Angle::from_arcsec)
        .def("to_deg", &astrometry::Angle::to_deg)
        .def("to_rad", &astrometry::Angle::to_rad);

    py::class_<astrometry::RightAscension, astrometry::Angle>(m, "RightAscension")
        .def(py::init<astrometry::Angle>());

    py::class_<astrometry::Declination, astrometry::Angle>(m, "Declination")
        .def(py::init<astrometry::Angle>());

    // --- Physics Module ---
    py::class_<physics::CartesianStateTyped<core::GCRF>>(m, "CartesianStateGCRF")
        .def_property_readonly("epoch", [](const physics::CartesianStateTyped<core::GCRF>& self) { return self.epoch; })
        .def_property_readonly("position", [](const physics::CartesianStateTyped<core::GCRF>& self) { return self.position.to_eigen_si(); })
        .def_property_readonly("velocity", [](const physics::CartesianStateTyped<core::GCRF>& self) { return self.velocity.to_eigen_si(); })
        .def("norm_au", [](const physics::CartesianStateTyped<core::GCRF>& self) { return self.position.norm().to_au(); });

    // --- Ephemeris Module ---
    py::enum_<astdyn::ephemeris::CelestialBody>(m, "CelestialBody")
        .value("SUN", astdyn::ephemeris::CelestialBody::SUN)
        .value("MERCURY", astdyn::ephemeris::CelestialBody::MERCURY)
        .value("VENUS", astdyn::ephemeris::CelestialBody::VENUS)
        .value("EARTH", astdyn::ephemeris::CelestialBody::EARTH)
        .value("MARS", astdyn::ephemeris::CelestialBody::MARS)
        .value("JUPITER", astdyn::ephemeris::CelestialBody::JUPITER)
        .value("SATURN", astdyn::ephemeris::CelestialBody::SATURN)
        .value("URANUS", astdyn::ephemeris::CelestialBody::URANUS)
        .value("NEPTUNE", astdyn::ephemeris::CelestialBody::NEPTUNE)
        .value("MOON", astdyn::ephemeris::CelestialBody::MOON)
        .export_values();

    py::class_<astdyn::ephemeris::PlanetaryEphemeris>(m, "PlanetaryEphemeris")
        .def(py::init<>()) // Add default constructor for Python to create instances
        .def("get_position", [](astdyn::ephemeris::PlanetaryEphemeris& self, astdyn::ephemeris::CelestialBody body, astdyn::time::EpochTDB t) {
            return self.getPosition(body, t).to_eigen_si(); // returns Eigen::Vector3d in meters
        })
        .def_static("set_provider", &astdyn::ephemeris::PlanetaryEphemeris::setGlobalProvider);

    // --- IO Module ---
    py::class_<io::MPCParser>(m, "MPCParser")
        .def_static("parse_line", &io::MPCParser::parse_line)
        .def_static("parse_file", &io::MPCParser::parse_file);

    py::class_<observations::OpticalObservation>(m, "OpticalObservation")
        .def(py::init<>())
        .def_readwrite("time", &observations::OpticalObservation::time)
        .def_readwrite("ra", &observations::OpticalObservation::ra)
        .def_readwrite("dec", &observations::OpticalObservation::dec)
        .def_readwrite("sigma_ra", &observations::OpticalObservation::sigma_ra)
        .def_readwrite("sigma_dec", &observations::OpticalObservation::sigma_dec)
        .def_readwrite("magnitude", &observations::OpticalObservation::magnitude);

    // --- OD Module ---
    py::class_<orbit_determination::GoodingIOD>(m, "GoodingIOD")
        .def(py::init<>())
        .def("compute", &orbit_determination::GoodingIOD::compute);

    py::class_<orbit_determination::GoodingIODResult>(m, "GoodingIODResult")
        .def_readwrite("success", &orbit_determination::GoodingIODResult::success)
        .def_readwrite("error_message", &orbit_determination::GoodingIODResult::error_message)
        .def_readwrite("solutions", &orbit_determination::GoodingIODResult::solutions);

    py::class_<orbit_determination::GoodingIODResult::Solution>(m, "GoodingSolution")
        .def_readwrite("state", &orbit_determination::GoodingIODResult::Solution::state)
        .def_readwrite("rho1", &orbit_determination::GoodingIODResult::Solution::rho1)
        .def_readwrite("rho3", &orbit_determination::GoodingIODResult::Solution::rho3);

    py::class_<orbit_determination::ExtendedKalmanFilter>(m, "EKF")
        .def(py::init<std::shared_ptr<propagation::Propagator>>())
        .def("update", &orbit_determination::ExtendedKalmanFilter::update);

    py::class_<orbit_determination::EKFResult>(m, "EKFResult")
        .def_readwrite("state", &orbit_determination::EKFResult::state)
        .def_readwrite("covariance", &orbit_determination::EKFResult::covariance)
        .def_readwrite("innovation", &orbit_determination::EKFResult::innovation);

    py::class_<orbit_determination::CovariancePropagator<core::GCRF>>(m, "CovariancePropagator")
        .def(py::init<std::shared_ptr<propagation::Propagator>>())
        .def("set_initial", &orbit_determination::CovariancePropagator<core::GCRF>::set_initial)
        .def("propagate", &orbit_determination::CovariancePropagator<core::GCRF>::propagate)
        .def("get_covariance", &orbit_determination::CovariancePropagator<core::GCRF>::get_covariance)
        .def("get_state", &orbit_determination::CovariancePropagator<core::GCRF>::get_state)
        .def("get_stm", &orbit_determination::CovariancePropagator<core::GCRF>::get_stm);

    // --- Propagation Module ---
    py::class_<propagation::PropagatorSettings>(m, "PropagatorSettings")
        .def(py::init<>())
        .def_readwrite("include_planets", &propagation::PropagatorSettings::include_planets)
        .def_readwrite("include_relativity", &propagation::PropagatorSettings::include_relativity)
        .def_readwrite("include_moon", &propagation::PropagatorSettings::include_moon)
        .def_readwrite("include_asteroids", &propagation::PropagatorSettings::include_asteroids)
        .def_readwrite("perturb_mercury", &propagation::PropagatorSettings::perturb_mercury)
        .def_readwrite("perturb_venus", &propagation::PropagatorSettings::perturb_venus)
        .def_readwrite("perturb_earth", &propagation::PropagatorSettings::perturb_earth)
        .def_readwrite("perturb_mars", &propagation::PropagatorSettings::perturb_mars)
        .def_readwrite("perturb_jupiter", &propagation::PropagatorSettings::perturb_jupiter)
        .def_readwrite("perturb_saturn", &propagation::PropagatorSettings::perturb_saturn)
        .def_readwrite("perturb_uranus", &propagation::PropagatorSettings::perturb_uranus)
        .def_readwrite("perturb_neptune", &propagation::PropagatorSettings::perturb_neptune)
        .def_readwrite("central_body_gm", &propagation::PropagatorSettings::central_body_gm)
        .def_readwrite("ppn_beta", &propagation::PropagatorSettings::ppn_beta)
        .def_readwrite("ppn_gamma", &propagation::PropagatorSettings::ppn_gamma)
        .def_readwrite("asteroid_ephemeris_file", &propagation::PropagatorSettings::asteroid_ephemeris_file)
        .def_readwrite("integrate_in_ecliptic", &propagation::PropagatorSettings::integrate_in_ecliptic)
        .def_readwrite("include_yarkovsky", &propagation::PropagatorSettings::include_yarkovsky)
        .def_readwrite("yarkovsky_a2", &propagation::PropagatorSettings::yarkovsky_a2);

    py::class_<propagation::Integrator, std::shared_ptr<propagation::Integrator>>(m, "Integrator")
        .def("reset_statistics", &propagation::Integrator::reset_statistics);

    py::class_<propagation::RK4Integrator, propagation::Integrator, std::shared_ptr<propagation::RK4Integrator>>(m, "RK4Integrator")
        .def(py::init<double>());

    py::class_<propagation::RKF78Integrator, propagation::Integrator, std::shared_ptr<propagation::RKF78Integrator>>(m, "RKF78Integrator")
        .def(py::init<double, double>());

    py::class_<propagation::AASIntegrator, propagation::Integrator, std::shared_ptr<propagation::AASIntegrator>>(m, "AASIntegrator")
        .def(py::init<double, double, double, double>(), 
             py::arg("precision") = 1e-4, py::arg("mu") = 1.32712440018e20, 
             py::arg("J2") = 0.0, py::arg("R_eq") = 6378137.0);

    py::class_<propagation::SABA4Integrator, propagation::Integrator, std::shared_ptr<propagation::SABA4Integrator>>(m, "SABA4Integrator")
        .def(py::init<double, double, double, double>(),
             py::arg("initial_step") = 1.0, py::arg("tolerance") = 1e-6,
             py::arg("min_step") = 1e-6, py::arg("max_step") = 100.0);

    py::class_<propagation::Propagator, std::shared_ptr<propagation::Propagator>>(m, "Propagator")
        .def(py::init<std::shared_ptr<propagation::Integrator>, std::shared_ptr<ephemeris::PlanetaryEphemeris>, const propagation::PropagatorSettings&>(),
             py::arg("integrator"), py::arg("ephem"), py::arg("settings") = propagation::PropagatorSettings())
        .def("propagate_cartesian", [](propagation::Propagator& self, const physics::CartesianStateTyped<core::GCRF>& initial, const time::EpochTDB& t) {
            return self.propagate_cartesian(initial, t);
        });

    // --- Ephemeris Providers ---
    py::class_<ephemeris::EphemerisProvider, std::shared_ptr<ephemeris::EphemerisProvider>>(m, "EphemerisProvider");
    
    py::class_<ephemeris::DE441Provider, ephemeris::EphemerisProvider, std::shared_ptr<ephemeris::DE441Provider>>(m, "DE441Provider")
        .def(py::init<const std::string&>());

    // --- High-Level APIs ---
    py::class_<AstDynEngine>(m, "AstDynEngine")
        .def(py::init<>())
        .def("load_observations", &AstDynEngine::load_observations)
        .def("fit_orbit", &AstDynEngine::fit_orbit)
        .def("compute_ephemeris", &AstDynEngine::compute_ephemeris);

    py::class_<api::OrbitFitAPI>(m, "OrbitFitAPI")
        .def_static("run_fit", &api::OrbitFitAPI::run_fit);

    m.def("analyze_orbit", [](const physics::CartesianStateTyped<core::GCRF>& state, 
                               const std::vector<observations::OpticalObservation>& obs) {
         // Default to non-perturbing if no provider set
         auto propagator = std::make_shared<propagation::Propagator>(
            std::make_unique<propagation::RKF78Integrator>(0.1),
            std::make_shared<ephemeris::PlanetaryEphemeris>(),
            propagation::PropagatorSettings{.include_planets = false}
         );
         auto summary = orbit_determination::ResidualAnalysis::analyze_orbit(state, obs, propagator);
         return summary.report_text;
    });

    // --- Occultation Module ---
    py::class_<astrometry::OccultationParameters>(m, "OccultationParameters")
        .def_readwrite("impact_parameter_km", &astrometry::OccultationParameters::impact_parameter_km)
        .def_readwrite("shadow_velocity_kms", &astrometry::OccultationParameters::shadow_velocity_kms)
        .def_readwrite("position_angle_deg", &astrometry::OccultationParameters::position_angle_deg)
        .def_readwrite("d_ra_cos_dec_mas_sec", &astrometry::OccultationParameters::d_ra_cos_dec_mas_sec)
        .def_readwrite("d_dec_mas_sec", &astrometry::OccultationParameters::d_dec_mas_sec)
        .def_readwrite("closest_approach_time_offset_sec", &astrometry::OccultationParameters::closest_approach_time_offset_sec);

    py::class_<astrometry::OccultationLogic>(m, "OccultationLogic")
        .def_static("compute_parameters", &astrometry::OccultationLogic::compute_parameters,
            py::arg("star_ra"), py::arg("star_dec"),
            py::arg("ast_ra"), py::arg("ast_dec"),
            py::arg("ast_dist_m"),
            py::arg("ast_dra_dt_rad_s"), py::arg("ast_ddec_dt_rad_s"),
            py::arg("ast_ddist_dt_m_s") = 0.0);

    m.attr("VERSION") = Version::string;
}
