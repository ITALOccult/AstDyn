#ifndef ASTDYN_PROPAGATOR_SETTINGS_HPP
#define ASTDYN_PROPAGATOR_SETTINGS_HPP

#include "astdyn/core/Constants.hpp"
#include "astdyn/core/IOCConfig.hpp"

#include <string>
#include <vector>

namespace astdyn::propagation {

struct PropagatorSettings {
    bool include_planets = true;
    bool include_moon = true;
    bool include_asteroids = true;
    std::vector<int> include_asteroids_list = {};
    std::vector<int> exclude_asteroids_list = {};
    bool use_default_asteroid_set = true;
    bool use_default_30_set = false;

    bool perturb_mercury = true;
    bool perturb_venus = true;
    bool perturb_earth = true;
    bool perturb_mars = true;
    bool perturb_jupiter = true;
    bool perturb_saturn = true;
    bool perturb_uranus = true;
    bool perturb_neptune = true;

    double central_body_gm = constants::GMS;

    bool include_relativity = true;
    double ppn_beta = 1.0;
    double ppn_gamma = 1.0;

    bool include_earth_j2 = true;
    bool include_sun_j2 = true;

    std::string asteroid_ephemeris_file = "";

    bool integrate_in_ecliptic = true;

    bool include_yarkovsky = false;
    double yarkovsky_a2 = 0.0;

    bool baricentric_integration = false;

    static PropagatorSettings from_config(const core::IOCConfig& cfg) {
        PropagatorSettings s;
        s.include_planets         = cfg.get<bool>  ("physics.planets",               s.include_planets);
        s.include_moon            = cfg.get<bool>  ("physics.moon",                  s.include_moon);
        s.include_asteroids       = cfg.get<bool>  ("physics.asteroids.enabled",     s.include_asteroids);
        s.use_default_30_set      = cfg.get<bool>  ("physics.asteroids.use_30_set",  s.use_default_30_set);
        s.perturb_mercury         = cfg.get<bool>  ("physics.mercury",               s.perturb_mercury);
        s.perturb_venus           = cfg.get<bool>  ("physics.venus",                 s.perturb_venus);
        s.perturb_earth           = cfg.get<bool>  ("physics.earth",                 s.perturb_earth);
        s.perturb_mars            = cfg.get<bool>  ("physics.mars",                  s.perturb_mars);
        s.perturb_jupiter         = cfg.get<bool>  ("physics.jupiter",               s.perturb_jupiter);
        s.perturb_saturn          = cfg.get<bool>  ("physics.saturn",                s.perturb_saturn);
        s.perturb_uranus          = cfg.get<bool>  ("physics.uranus",                s.perturb_uranus);
        s.perturb_neptune         = cfg.get<bool>  ("physics.neptune",               s.perturb_neptune);
        s.include_relativity      = cfg.get<bool>  ("physics.relativity",            s.include_relativity);
        s.ppn_beta                = cfg.get<double>("physics.ppn_beta",              s.ppn_beta);
        s.ppn_gamma               = cfg.get<double>("physics.ppn_gamma",             s.ppn_gamma);
        s.include_earth_j2        = cfg.get<bool>  ("physics.earth_j2",             s.include_earth_j2);
        s.include_sun_j2          = cfg.get<bool>  ("physics.sun_j2",               s.include_sun_j2);
        s.include_yarkovsky       = cfg.get<bool>  ("physics.yarkovsky.enabled",     s.include_yarkovsky);
        s.yarkovsky_a2            = cfg.get<double>("physics.yarkovsky.a2",          s.yarkovsky_a2);
        s.baricentric_integration = cfg.get<bool>  ("integrator.barycentric",        s.baricentric_integration);
        s.integrate_in_ecliptic   = cfg.get<bool>  ("integrator.ecliptic_frame",     s.integrate_in_ecliptic);
        return s;
    }
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATOR_SETTINGS_HPP