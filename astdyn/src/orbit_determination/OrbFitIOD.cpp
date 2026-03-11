#include "astdyn/orbit_determination/OrbFitIOD.hpp"
#include "astdyn/orbfit/iod.h"
#include "astdyn/orbfit/observations.h"
#include "astdyn/math/frame_algebra.hpp"
#include "src/core/frame_tags.hpp"
#include <iostream>
#include <algorithm>

#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/observations/ObservatoryDatabase.hpp"

namespace astdyn::orbit_determination {

OrbFitIOD::OrbFitIOD(const GaussIODSettings& settings) : settings_(settings) {}

std::optional<std::array<int, 3>> OrbFitIOD::select_observations(
    const std::vector<observations::OpticalObservation>& observations) const {
    if (observations.size() < 3) return std::nullopt;

    int n = observations.size();
    return std::array<int, 3>{0, n / 2, n - 1};
}

GaussIODResult OrbFitIOD::compute_from_three(
    const observations::OpticalObservation& obs1,
    const observations::OpticalObservation& obs2,
    const observations::OpticalObservation& obs3) {

    auto convert_obs = [](const observations::OpticalObservation& o) -> orbfit::Observation {
        orbfit::Observation out;
        out.t_obs = o.time.mjd() + 2400000.5; // JD
        out.ra = o.ra.to_rad();
        out.dec = o.dec.to_rad();
        out.sigma_ra = o.sigma_ra.to_rad();
        out.sigma_dec = o.sigma_dec.to_rad();
        out.station_code = o.observatory_code;
        
        auto t_tdb = astdyn::time::EpochTDB::from_mjd(o.time.mjd()); // Just approximate UTC=TDB for this
        auto earth = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t_tdb);
        math::Vector3<core::GCRF, physics::Distance> R = earth.position;
        const auto& obs_db = observations::ObservatoryDatabase::getInstance();
        if (auto info = obs_db.getObservatory(o.observatory_code)) {
            auto topo = info->getPositionGCRF(o.time);
            R = R + math::Vector3<core::GCRF, physics::Distance>::from_si(topo.x_si(), topo.y_si(), topo.z_si());
        }
        
        auto dummy_state = physics::CartesianStateTyped<core::GCRF>::from_si(
            t_tdb, R.x_si(), R.y_si(), R.z_si(), 0.0, 0.0, 0.0
        );
        auto state_eclip = dummy_state.cast_frame<core::ECLIPJ2000>();
        auto r_eclip = state_eclip.position;
        out.obs_pos = orbfit::Vec3(r_eclip.x_si() / constants::AU, r_eclip.y_si() / constants::AU, r_eclip.z_si() / constants::AU);
        
        return out;
    };

    orbfit::Observation o1 = convert_obs(obs1);
    orbfit::Observation o2 = convert_obs(obs2);
    orbfit::Observation o3 = convert_obs(obs3);

    GaussIODResult result;
    try {
        orbfit::GaussSolution sol;
        double dt13 = std::abs(o3.t_obs - o1.t_obs);
        if (dt13 < 2.0) {
            sol = orbfit::laplace_iod(o1, o2, o3);
        } else {
            sol = orbfit::gauss_iod(o1, o2, o3);
        }

        result.success = sol.converged;
        result.error_message = "";
        
        // OrbFit returns state in heliocentric ECLIPJ2000
        orbfit::StateVector sv = orbfit::elements_to_state(sol.elements, sol.elements.epoch);
        
        auto epoch = time::EpochTDB::from_mjd(sol.elements.epoch - 2400000.5);
        Eigen::VectorXd y(6);
        y << sv.pos.x, sv.pos.y, sv.pos.z, sv.vel.x, sv.vel.y, sv.vel.z;
        auto state_eclip = physics::CartesianStateTyped<core::ECLIPJ2000>::from_au_aud(
            epoch, y
        );

        result.state = state_eclip.cast_frame<core::GCRF>();
        result.epoch = epoch;
        
        // Dummy slant ranges as OrbFit's API doesn't expose them directly in the return struct
        result.slant_range_1 = (sv.pos - o1.obs_pos).norm();
        result.slant_range_2 = (sv.pos - o2.obs_pos).norm();
        result.slant_range_3 = (sv.pos - o3.obs_pos).norm();
        
        result.iterations = 1;

    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = e.what();
    }

    return result;
}

GaussIODResult OrbFitIOD::compute(const std::vector<observations::OpticalObservation>& observations) {
    auto maybe_indices = select_observations(observations);
    if (!maybe_indices) {
        GaussIODResult res;
        res.success = false;
        res.error_message = "Not enough observations";
        return res;
    }

    auto indices = *maybe_indices;
    auto res = compute_from_three(observations[indices[0]], observations[indices[1]], observations[indices[2]]);
    res.obs_index_1 = indices[0];
    res.obs_index_2 = indices[1];
    res.obs_index_3 = indices[2];
    return res;
}

} // namespace astdyn::orbit_determination
