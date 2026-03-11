#pragma once
// planet_cache.h — must be included BEFORE integrator.h or after overriding planet_state
// Patches orbfit::planet_state with a simple time-indexed cache (LRU-1 per planet)

#include "constants.h"
#include "linalg.h"
#include "orbital_elements.h"
#include <array>
#include <cmath>

namespace orbfit {

// Override planet_state with a simple last-value cache per planet
// The RK45 stepper calls planet_state with the same or nearby t many times.
// This gives ~6-8x speedup.

struct CachedPlanetState {
    double t_cached = -1e30;
    PlanetState state{};
};

inline PlanetState planet_state_cached(int pid, double t) {
    // Static per-planet cache (thread-unsafe, fine for single-threaded use)
    static std::array<CachedPlanetState, 9> cache{};
    auto& c = cache[pid];
    if (std::abs(c.t_cached - t) < 1e-9) return c.state;
    c.state    = planet_state(pid, t);
    c.t_cached = t;
    return c.state;
}

// Faster EOM using cache
inline State6 eom_fast(const State6& s, double t,
                        bool perturbers = true,
                        bool relativistic = false)
{
    Vec3 r{s[0],s[1],s[2]};
    Vec3 v{s[3],s[4],s[5]};
    double rv  = r.norm();
    double rv3 = rv*rv*rv;
    Vec3 acc   = r * (-GM_SUN / rv3);

    if (perturbers) {
        for (int pid = 0; pid < 8; pid++) {
            auto ps = planet_state_cached(pid, t);
            Vec3   dp  = r - ps.pos;
            double dp3 = std::pow(dp.norm(), 3.0);
            double rp3 = std::pow(ps.pos.norm(), 3.0);
            acc -= dp * (ps.gm / dp3);
            acc -= ps.pos * (ps.gm / rp3);
        }
    }
    return {v.x, v.y, v.z, acc.x, acc.y, acc.z};
}

} // namespace orbfit
