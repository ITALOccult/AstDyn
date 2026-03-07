#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "astdyn/AstDyn.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/io/MPCParser.hpp"
#include "astdyn/orbit_determination/GoodingIOD.hpp"
#include "astdyn/orbit_determination/ResidualAnalysis.hpp"

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

    // --- Astrometry Module ---
    py::class_<astrometry::Angle>(m, "Angle")
        .def_static("from_deg", &astrometry::Angle::from_deg)
        .def_static("from_rad", &astrometry::Angle::from_rad)
        .def_static("from_arcsec", &astrometry::Angle::from_arcsec)
        .def("to_deg", &astrometry::Angle::to_deg)
        .def("to_rad", &astrometry::Angle::to_rad);

    // --- Physics Module ---
    py::class_<physics::CartesianStateTyped<core::GCRF>>(m, "CartesianStateGCRF")
        .def_property_readonly("epoch", [](const physics::CartesianStateTyped<core::GCRF>& self) { return self.epoch; })
        .def_property_readonly("position", [](const physics::CartesianStateTyped<core::GCRF>& self) { return self.position.to_eigen_si(); })
        .def_property_readonly("velocity", [](const physics::CartesianStateTyped<core::GCRF>& self) { return self.velocity.to_eigen_si(); })
        .def("norm_au", [](const physics::CartesianStateTyped<core::GCRF>& self) { return self.position.norm().to_au(); });

    // --- IO Module ---
    py::class_<io::MPCParser>(m, "MPCParser")
        .def_static("parse_line", &io::MPCParser::parse_line)
        .def_static("parse_file", &io::MPCParser::parse_file);

    py::class_<observations::OpticalObservation>(m, "OpticalObservation")
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

    m.def("analyze_orbit", [](const physics::CartesianStateTyped<core::GCRF>& state, 
                              const std::vector<observations::OpticalObservation>& obs) {
         auto propagator = std::make_shared<propagation::Propagator>(
            std::make_unique<propagation::RKF78Integrator>(0.1),
            std::make_shared<ephemeris::PlanetaryEphemeris>()
         );
         auto summary = orbit_determination::ResidualAnalysis::analyze_orbit(state, obs, propagator);
         return summary.report_text;
    });
}
